import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Created by Геннадий on 08.12.2014.
 */
class ProblemIntervalBuilder extends Constants{

  private Reference genome;
  private int windowSize;
  private MappedData mappedData;
  private float identityThreshold;
  private int K;
  private int anchorSize;
  private Logger logger;
  private int threadNum;

  ProblemIntervalBuilder(Reference genome, MappedData mappedData,
                                int windowSize, float identityThreshold, int K,
                                int anchorSize, Logger logger, int threadNum){
    this.genome = genome;
    this.windowSize = windowSize;
    this.mappedData = mappedData;
    this.identityThreshold = identityThreshold;
    this.K = K;
    this.anchorSize = anchorSize;
    this.logger = logger;
    this.threadNum = threadNum;
  }

  private class IdentityCounting implements Runnable{

    public float[] sumIdentity;
    public int[] readCount;
    public int from;
    public int to;


    @Override
    public void run(){
      sumIdentity = new float[genome.length];
      readCount = new int[genome.length];
      Iterator<MappedRead> iter = mappedData.mappedReads.listIterator(from);
      for(int n = from; n < to; n++){
        MappedRead read = iter.next();
        if(read.seq.length < windowSize){
          continue;
        }
        Alignment alignment = read.aln;
        byte[] seq1 = alignment.sequence1;
        int identity = 0;
        int startNoClip = read.start + alignment.start2;
        int start = Math.max(0, startNoClip - alignment.start1);
        int endNoClip = read.start + alignment.end2;
        int end = Math.min(genome.length, endNoClip + read.seq.length - alignment.end1);
        int alnEndPos = 0;
        int len;
        int alnStartPos = 0;
        if(alignment.start1 >= windowSize){
          len = windowSize;
        }else{
          len = alignment.start1;
          if(alignment.length > windowSize - alignment.start1){
            for(int i = startNoClip; i < startNoClip + windowSize - alignment.start1 &&
                alnEndPos < alignment.length;){
              byte c = seq1[alnEndPos];
              if(c != GAP && c == genome.seqB[i]){
                identity ++;
                alnEndPos++;
                len ++;
                i++;
              }else{
                if(c < 'a'){
                  alnEndPos++;
                  len ++;
                  i++;
                }else{
                  do{
                    alnEndPos ++;
                    len ++;
                  }while(seq1[alnEndPos] >= 'a');
                }
              }
            }
            while(alnEndPos < alignment.length && seq1[alnEndPos] >= 'a'){
              alnEndPos ++;
              len ++;
            }
          }else{
            identity += alignment.identity;
            len += windowSize - alignment.start1;
            alnEndPos = alignment.length;
          }
        }
        sumIdentity[start] += (float)identity/len;
        readCount[start] ++;
        for(int i = start + 1; i <= end - windowSize; i++){
          if(i > startNoClip && i <= endNoClip){
            byte c = seq1[alnStartPos];
            if(c != GAP && c == genome.seqB[i - 1]){
              identity --;
            }
            alnStartPos ++;
            if(alnStartPos < alignment.length){
              while(seq1[alnStartPos] >= 'a'){
                alnStartPos ++;
                len --;
              }
            }
          }
          if(i > startNoClip - windowSize && i <= endNoClip - windowSize){
            byte c = seq1[alnEndPos];
            if(c != GAP && c == genome.seqB[i + windowSize - 1]){
              identity ++;
            }
            alnEndPos ++;
            if(alnEndPos < alignment.length){
              while(seq1[alnEndPos] >= 'a'){
                alnEndPos ++;
                len ++;
              }
            }
          }
          sumIdentity[i] += (float)identity/len;
          readCount[i] ++;
        }
      }
    }
  }

  private void countIdentity(float[] identity, int[] readCount)
      throws InterruptedException{
    IdentityCounting[] tasks = new IdentityCounting[threadNum];
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    int size = mappedData.mappedReads.size()/threadNum;
    for(int i = 0; i < threadNum - 1; i++){
      IdentityCounting task = new IdentityCounting();
      task.from = i*size;
      task.to = task.from + size;
      executor.execute(task);
      tasks[i] = task;
    }
    IdentityCounting task = new IdentityCounting();
    task.from = (threadNum - 1)*size;
    task.to = mappedData.mappedReads.size();
    executor.execute(task);
    tasks[threadNum - 1] = task;
    executor.shutdown();
    executor.awaitTermination(10, TimeUnit.HOURS);
    for(IdentityCounting t: tasks){
      for(int i = 0; i < genome.length; i++){
        identity[i] += t.sumIdentity[i];
        readCount[i] += t.readCount[i];
      }
    }
  }

  Interval[] localizeProblemRegions()
      throws IOException, InterruptedException{
    float[] identity = new float[genome.length];
    int[] readCounts = new int[genome.length];
    countIdentity(identity, readCounts);
    boolean isGap = false;
    ArrayList<Interval> problemRegions = new ArrayList<>();
    Interval curr = new Interval();
    for(int i = 0; i < genome.length; i++){
      if(genome.seqB[i] == GAP){
        if(!isGap){
          isGap = true;
          curr.start = i;
        }
      }else{
        if(isGap){
          isGap = false;
          curr.end = i;
          problemRegions.add(curr);
          curr = new Interval();
        }
        float avIdentity;
        if(readCounts[i] > 0){
          avIdentity = identity[i]/readCounts[i];
        }else{
          avIdentity = 0;
        }
        if(avIdentity < identityThreshold){
          Interval interval = new Interval();
          interval.start = i;
          interval.end = i + windowSize;
          problemRegions.add(interval);
        }
      }
    }
    if(isGap){
      curr.end = genome.length;
      problemRegions.add(curr);
    }
    ArrayList<Interval> problemIntervalsArr = new ArrayList<>();
    if(problemRegions.size() > 0){
      Collections.sort(problemRegions);
      Iterator<Interval> iter = problemRegions.iterator();
      Interval prev = iter.next();
      while(iter.hasNext()){
        curr = iter.next();
        if(curr.start - prev.end < anchorSize){
          prev.end = Math.max(prev.end, curr.end);
        }else{
          problemIntervalsArr.add(prev);
          prev = curr;
        }
      }
      if(prev.end > genome.length){
        prev.end = genome.length;
      }
      problemIntervalsArr.add(prev);
    }

    ArrayList<Interval> res = new ArrayList<>();
    if(genome.isFragmented){
      for(Interval interval: problemIntervalsArr){
        for(int i = 0; i < genome.fragmentEnds.length; i++){
          int end_i = genome.fragmentEnds[i];
          if(interval.start < end_i){
            if(interval.end > end_i){
              Interval leftInterval = new Interval();
              leftInterval.start = interval.start;
              leftInterval.end = end_i;
              res.add(leftInterval);
              int endPrev = end_i;
              for(int j = i + 1; j < genome.fragmentEnds.length; j++){
                int end_j = genome.fragmentEnds[j];
                if(interval.end <= end_j){
                  interval.start = endPrev;
                  break;
                }else{
                  Interval middleInterval = new Interval();
                  middleInterval.start = endPrev;
                  middleInterval.end = end_j;
                  res.add(middleInterval);
                  endPrev = end_j;
                }
              }
            }
            break;
          }
        }
        res.add(interval);
      }
      for(Interval interval: res){
        int fragmentStart = 0;
        for(int i = 0; i < genome.fragmentEnds.length; i++){
          int fragmentEnd = genome.fragmentEnds[i];
          if(interval.end <= fragmentEnd){
            if(interval.start - fragmentStart < anchorSize){
              if(interval.end > fragmentEnd - anchorSize){
                interval.seqFinal = "------";
              }else{
                interval.start = fragmentStart;
                interval.concat = new byte[anchorSize];
                interval.left = new byte[0];
                interval.right = new byte[anchorSize];
                for(int j = 0, n = interval.end; j < anchorSize; j++, n++){
                  interval.concat[j] = genome.seqB[n];
                  interval.right[j] = genome.seqB[n];
                }
                interval.junction = -1;
                interval.init(K, threadNum);
              }
            }else{
              if(interval.end > fragmentEnd - anchorSize){
                interval.concat = new byte[anchorSize];
                interval.left = new byte[anchorSize];
                interval.right = new byte[0];
                interval.end = fragmentEnd;
                for(int j = 0, n = interval.start - anchorSize; j < anchorSize; j++, n++){
                  interval.concat[j] = genome.seqB[n];
                  interval.left[j] = genome.seqB[n];
                }
                interval.junction = (short) (anchorSize - 1);
                interval.init(K, threadNum);
              }else{
                interval.concat = new byte[2*anchorSize];
                interval.left = new byte[anchorSize];
                interval.right = new byte[anchorSize];
                for(int m = 0, j = interval.start - anchorSize, k = anchorSize, n = interval.end;
                    m < anchorSize; m++, j++, k++, n++){
                  interval.concat[m] = genome.seqB[j];
                  interval.left[m] = genome.seqB[j];
                  interval.concat[k] = genome.seqB[n];
                  interval.right[m] = genome.seqB[n];
                }
                interval.junction = (short) (anchorSize - 1);
                interval.init(K, threadNum);
              }
            }
            break;
          }
          fragmentStart = fragmentEnd;
        }
      }
    }else{
      res = problemIntervalsArr;
      for(Interval interval: problemIntervalsArr){
        if(interval.start < anchorSize){
          if(interval.end > genome.length - anchorSize){
            logger.printf("Not enough reads to assemble\n");
            interval.seqFinal = "------";
          }else{
            interval.start = 0;
            interval.concat = new byte[anchorSize];
            interval.left = new byte[0];
            interval.right = new byte[anchorSize];
            for(int i = 0, n = interval.end; i < anchorSize; i++, n++){
              interval.concat[i] = genome.seqB[n];
              interval.right[i] = genome.seqB[n];
            }
            interval.junction = -1;
            interval.init(K, threadNum);
          }
        }else{
          if(interval.end > genome.length - anchorSize){
            interval.concat = new byte[anchorSize];
            interval.left = new byte[anchorSize];
            interval.right = new byte[0];
            interval.end = genome.length;
            for(int i = 0, n = interval.start - anchorSize; i < anchorSize; i++, n++){
              interval.concat[i] = genome.seqB[n];
              interval.left[i] = genome.seqB[n];
            }
            interval.junction = (short) (anchorSize - 1);
            interval.init(K, threadNum);
          }else{
            interval.concat = new byte[2*anchorSize];
            interval.left = new byte[anchorSize];
            interval.right = new byte[anchorSize];
            for(int i = 0, j = interval.start - anchorSize, k = anchorSize, n = interval.end;
                i < anchorSize; i++, j++, k++, n++){
              interval.concat[i] = genome.seqB[j];
              interval.left[i] = genome.seqB[j];
              interval.concat[k] = genome.seqB[n];
              interval.right[i] = genome.seqB[n];
            }
            interval.junction = (short) (anchorSize - 1);
            interval.init(K, threadNum);
          }
        }
      }
    }
    return res.toArray(new Interval[res.size()]);
  }


}
