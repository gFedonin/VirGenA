import org.jdom2.Document;
import org.jdom2.Element;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Created by Геннадий on 09.12.2014.
 */
class ConsensusBuilderSimple extends ConsensusBuilder{

  MappedData mappedData;
  Mapper mapper;
  int coverageThreshold;
  int threadNum;

  ConsensusBuilderSimple(Document document){
    Element element = document.getRootElement().getChild("ConsensusBuilder");
    coverageThreshold = Integer.parseInt(element.getChildText("CoverageThreshold"));
    mapper = new Mapper(document);
    threadNum = Integer.parseInt(document.getRootElement().getChildText("ThreadNumber"));
    if(threadNum == -1){
      threadNum = Runtime.getRuntime().availableProcessors();
    }
  }

  @Override
  public String buildConsensus(Reference genome, ArrayList<PairedRead> reads){
    try{
      mappedData = mapper.mapReads(reads, genome);
      return new String(getConsensusSimple(false, coverageThreshold, genome.seqB));
    }catch(Exception e){
      e.printStackTrace();
    }
    return null;
  }

  private class ReadCounting implements Runnable{

    public int[][] counts;
    public int from;
    public int to;
    public int genomeLen;
    //private boolean concordantOnly;

    public ReadCounting(int from, int to, int genomeLen){//, boolean concordantOnly
        this.from = from;
        this.to = to;
        this.genomeLen = genomeLen;
        //this.concordantOnly = concordantOnly;
    }

    @Override
    public void run(){
      counts = new int[genomeLen][5];
      Iterator<MappedRead> iter = mappedData.mappedReads.listIterator(from);
      for(int n = from; n < to; n++){
        MappedRead read = iter.next();
//        if(mappedData.isConcordant[n] == 0){
//          continue;
//        }
        Alignment alignment = read.aln;
        byte[] seq1 = alignment.sequence1;
        int j = read.start + alignment.start2;
        for(int i = 0; i < alignment.length; i++){
          if(seq1[i] < (byte)'a'){
            counts[j][nToI[seq1[i]]] ++;
            j ++;
          }
        }
      }
    }
  }

  private class ConsensusBuildingSimple implements Runnable{

    int[][] counts;
    byte[] consensus;
    private int from;
    private int to;
    private int coverageThreshold;
    private boolean saveUncovered;
    private byte[] genome;

    ConsensusBuildingSimple(int from, int to, int coverageThreshold,
                                   boolean saveUncovered, byte[] genome, int[][] counts){
        this.from = from;
        this.to = to;
        this.coverageThreshold = coverageThreshold;
        this.saveUncovered = saveUncovered;
        this.genome = genome;
        this.counts = counts;
    }

    @Override
    public void run(){
      consensus = new byte[to - from];
      for(int i = from; i < to; i++){
        int maxCount = counts[i][0];
        int maxCountIndex = 0;
        for(int j = 1; j < 4; j++){
          if(counts[i][j] > maxCount){
            maxCount = counts[i][j];
            maxCountIndex = j;
          }
        }
        if(maxCount <= coverageThreshold){
          if(saveUncovered){
            consensus[i - from] = genome[i];
          }else{
            consensus[i - from] = GAP;
          }
        }else{
          consensus[i - from] = iToN[maxCountIndex];
        }
      }
    }
  }

  byte[] getConsensusSimple(boolean saveUncovered, int coverageThreshold, byte[] genome) throws InterruptedException{//boolean concordantOnly,
    ReadCounting[] tasks = new ReadCounting[threadNum];
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    int size = mappedData.mappedReads.size()/threadNum;
    for(int i = 0; i < threadNum - 1; i++){
      ReadCounting task = new ReadCounting(i*size, (i + 1)*size, genome.length);
      executor.execute(task);
      tasks[i] = task;
    }
    ReadCounting task = new ReadCounting((threadNum - 1)*size, mappedData.mappedReads.size(),
            genome.length);
    executor.execute(task);
    tasks[threadNum - 1] = task;
    executor.shutdown();
    executor.awaitTermination(10, TimeUnit.HOURS);
    int[][] counts = new int[genome.length][5];
    for(ReadCounting t: tasks){
      for(int i = 0; i < genome.length; i++){
        for(int j = 0; j < 5; j++){
          counts[i][j] += t.counts[i][j];
        }
      }
    }
    tasks = null;
    byte[] consensus = new byte[genome.length];
    ConsensusBuildingSimple[] tasksB = new ConsensusBuildingSimple[threadNum];
    executor = Executors.newFixedThreadPool(threadNum);
    size = genome.length/threadNum;
    for(int i = 0; i < threadNum - 1; i++){
      ConsensusBuildingSimple taskB = new ConsensusBuildingSimple(i*size, (i + 1)*size,
              coverageThreshold, saveUncovered, genome, counts);
      executor.execute(taskB);
      tasksB[i] = taskB;
    }
    ConsensusBuildingSimple taskB = new ConsensusBuildingSimple((threadNum - 1)*size,
            genome.length, coverageThreshold, saveUncovered, genome, counts);
    executor.execute(taskB);
    tasksB[threadNum - 1] = taskB;
    executor.shutdown();
    executor.awaitTermination(10, TimeUnit.HOURS);
    for(ConsensusBuildingSimple t: tasksB){
      System.arraycopy(t.consensus, 0, consensus, t.from, t.to - t.from);
    }
    return consensus;
  }


}
