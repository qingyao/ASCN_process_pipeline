##aroma-final##
import click, sys, os, time, logging, socket, subprocess
import pandas as pd
from subprocess import run
from multiprocessing.dummy import Pool
import finder_colors as color

@click.command()
@click.option('-b','--block', default=0, help='Block number out of 5 blocks: 0--4')
@click.option('-c','--cleanup', default=1, help='Clean up intermediate files: 1 for yes, 0 for no')
@click.option('-t','--threads', default=4, help='Number of threads to use')
@click.option('-w','--workingdir', default='/Users/bgprocess/aroma/hg19/', help='Aroma working directory')
@click.option('-s','--sourcedir',default='/Users/bgprocess/aroma/AromaPack/',help="Path to folder 'Aromapack'")
@click.option('-a','--allblocks',default=5,help='Number of machines to pass to')
@click.option('-m','--memory',default=100,help='memory setup for aroma pipelines, for <= 64GB machine use 50')
@click.option('-e','--error',default=0,help='if 1, will run the series which contains error')
@click.option('-f','--force',default=0,help='if 1, will force rerun for selected series')
@click.option('-o','--option',default='noprocess',help="'probe' for raw data to probes; 'seg' for initial segmentation; 'reseg' for evaluation of cn segments and possible re-segmentation")
@click.option('-ft','--filetype',default='cn',help="in probe processing not needed, for seg or reseg, need to indicate 'cn' or 'fracb'.")
def cli(block,cleanup,threads,workingdir,sourcedir,allblocks,memory,error,force,option,filetype):
    if option == 'noprocess':
        sys.exit('no process option specified.')
    p = Pool(processes=threads)
    serieslist=[]
    if error == 1:
        dir = os.path.join(workingdir,'processed')
        files = os.listdir(dir)
        files = [f for f in files if f.startswith('aroma')]
        print(files)
        st = set()
        for f in files:
            if os.stat(os.path.join(dir,f)).st_size != 0:
                aroma = pd.read_table(os.path.join(dir,f),header=None)

                e = aroma.iloc[:,3]
                for i in e:
                    if type(i).__name__ == 'str':
                        st.add(i)
        with open(os.path.join(sourcedir,'affy.tsv'),'w',encoding='utf-8') as wf:
            wf.write('\n'.join(st))
            wf.close()

    with open(os.path.join(sourcedir,'affy.tsv'),'r',encoding='utf-8') as f:
        c=0
        for line in f:
            if c%allblocks==block:
                if not line.startswith('G'):
                    serieslist.append('GSE'+line.rstrip())
                else:
                    serieslist.append(line.rstrip())
            c+=1
    # if not force:
    #     rmlist = []
    #     for s in serieslist:
    #         if os.path.exists('/Volumes/arraymapMirror/arraymap/hg19/{0}'.format(s)):
    #             if color.get('/Volumes/arraymapMirror/arraymap/hg19/{0}'.format(s)) == 'gray':
    #                 rmlist.append(s)
    #     serieslist = list(set(serieslist) - set(rmlist))

    machine = socket.gethostname()
    startTime = time.time()
    logfilename = '%s/%s/%sH.log' % (workingdir,'processed',time.strftime('%Y-%m-%d,%H'))
    formatter = logging.Formatter('%(asctime)s\t%(message)s')
    logger = logging.getLogger('process')
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(logfilename)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    logger.info ('Doing %d series: \n%s' %(len(serieslist),serieslist))

    def runPipeline(seriesname):
        seriespath = os.path.join(workingdir,'processed',seriesname)
        if not os.path.exists(seriespath):
            os.makedirs(seriespath)

        logger.debug('Machine: %s\tSeries: %s'%(machine, seriesname))
        logger.info('Started series: '+seriesname)
        seriesStart = time.time()
        if option == 'probe':
            shellcommand = 'R --vanilla <{0}/ProcessPipeline-probe.R --args {1} {2} {3} {0} {4} {5} &>/dev/null'.format(sourcedir,workingdir,seriesname,cleanup,memory,force)
        elif option == 'seg':
            shellcommand = 'R --vanilla <{0}/test_dir/ProcessPipeline-seg_test.R --args {1} {2} {0} {3} {4} &>/dev/null'.format(sourcedir,workingdir,seriesname,force,filetype)
        elif option == 'reseg':
            shellcommand = 'R --vanilla <{0}/ProcessPipeline-reseg.R --args {1} {2} {0} {3} &>/dev/null'.format(sourcedir,workingdir,seriesname, filetype)

        run(shellcommand,shell=True)

        seriesTime =round(time.time()-seriesStart)
        elapsedTime =round(time.time()-startTime)
        sm, ss = divmod(seriesTime, 60)
        sh, sm = divmod(sm, 60)
        em, es = divmod(elapsedTime, 60)
        eh, em = divmod(em, 60)
        logger.info('Finished series:'+seriesname+'. Used:'+"%dhr%02dmin%02dsec" % (sh, sm, ss)+'. Total elapsed:'+"%dhr%02dmin%02dsec" % (eh, em, es)+'.\n')
        # with open(os.path.join(seriespath,'logfile'),'a+',encoding='utf-8') as wf:
        #     wf.write('\t'.join([time.strftime('%Y-%m-%d %H:%M:%S'),machine,str(elapsedTime)+'sec\n']))
        #     wf.close()

        syncCmd = 'expect {0}/test_dir/expect_syncback_test.sh {1}'.format(sourcedir,seriesname,workingdir)
        rmCmd = 'for dir in {0}/processed/{1}/GSM*;do rm $dir/fracB,*;rm $dir/probes,*;done &>/dev/null'.format(workingdir,seriesname)
        run(syncCmd,shell=True)

        try:
            if not checksum(sourcedir,seriesname):
                logger.info('series {0} sync not complete. try again...'.format(seriesname))
                run(syncCmd,shell=True)
                if checksum(sourcedir,seriesname):
                    logger.info('series {0} sync complete in second try'.format(seriesname))
                    run(rmCmd,shell=True)
                else: logger.info('series {0} sync not complete after 2 trials'.format(seriesname))
            else:
                run(rmCmd,shell=True)
        except Exception as err:
            logger.info('dryrun error in series {0}: {1}: {2}'.format(seriesname,type(err),err))

    p.map(runPipeline,serieslist)

    shellcommand = 'python3 <{0}/set_color.py &>/dev/null'.format(sourcedir)
    run(shellcommand,shell=True)

def local_exists(path,searched):
    list=[]
    for _, _, files in os.walk(path):
        list+=files
    if searched in list:
        return True
    else:
        return False

def checksum(sourcedir,series):
    cmd_expectR = 'expect ' + sourcedir + '/dryrun.sh ' + series
    output = subprocess.check_output(cmd_expectR,shell=True)
    output = output.decode('utf-8')
    output = output.splitlines()
    if len(output) > 15: return False
    stidx = None
    for i,p in enumerate(output):
        if p=='building file list ... done' or p=='done' or p == 'sending incremental file list':
            stidx = i
            break
    if not stidx:
        print(output)
        return False
    elif output[stidx+1] == '' and output[stidx+2].startswith('sent'):
        return True
    elif output[stidx+1].startswith('created') and output[stidx+2] == '' and output[stidx+3].startswith('sent'):
        return True
    elif output[stidx+2].startswith('created') and output[stidx+3] == '' and output[stidx+4].startswith('sent'):
        return True
    else:
        return False

if __name__ == '__main__':
    print()
    cli()
