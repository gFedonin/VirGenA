import org.jdom2.Document;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Created by Gennady on 07.11.2015.
 */
class MapperMSA extends Constants{

  private int insertLen;
  private KMerCounterMSA counterRF;
  private int threadNum;
  Logger logger;
  private int batchSize;

  MapperMSA(Document document){
    counterRF = KMerCounterMSA.getInstance(document);
    logger = Logger.getInstance(document);
    insertLen = Integer.parseInt(document.getRootElement().getChild("Data").getChildText("InsertionLength"));
    threadNum = Integer.parseInt(document.getRootElement().getChildText("ThreadNumber"));
    batchSize = Integer.parseInt(document.getRootElement().getChildText("BatchSize"));
  }

  private class ReadMapping implements Runnable{

    public ArrayList<PairedRead> reads;
    public ReferenceAlignment refAlignment;
    public int from;
    public int to;
    public ArrayList<MappedRead> mappedReads;
    public ArrayList<PairedRead> concordant;
    public ArrayList<PairedRead> discordant;
    public ArrayList<PairedRead> leftMateMapped;
    public ArrayList<PairedRead> rightMateMapped;
    public ArrayList<PairedRead> unmapped;
    public int readsTotal = 0;
    public int bothExist = 0;
    public int mappedForward = 0;
    public int mappedReverse = 0;
    public long totalScore = 0;

    public ReadMapping(ArrayList<PairedRead> reads, ReferenceAlignment refAlignment,
                       int from, int to){
      this.reads =reads;
      this.refAlignment = refAlignment;
      this.from = from;
      this.to = to;
    }

    private MappedRead[] mapReadFast(byte[] read){
      if(Arrays.equals(read, NULL_SEQ)){
        return new MappedRead[0];
      }
      int readLen = read.length;
      MappedRead[] bestPositions = counterRF.getNBestRegions(read, refAlignment.refAlnIndex, 0,
          refAlignment.length, counterRF.cutoffs[readLen]);
      if(bestPositions == null){
        return new MappedRead[0];
      }
      for(MappedRead bestPosition : bestPositions){
        bestPosition.seq = read;
      }
      return bestPositions;
    }

    private void mapRead(PairedRead read){
      MappedRead[] forward1 = mapReadFast(read.seq1);
      MappedRead[] reverse1 = mapReadFast(getComplement(read.seq1));
      MappedRead[] forward2 = mapReadFast(read.seq2);
      MappedRead[] reverse2 = mapReadFast(getComplement(read.seq2));
      int maxCount = 0;
      int max1 = -1;
      int max2 = -1;
      byte side = -1;
      int maxInd1 = -1;
      int count1 = 0;
      int side1 = -1;
      int maxInd2 = -1;
      int count2 = 0;
      int side2 = -1;
      for(int j = 0; j < forward1.length; j++){
        MappedRead read1 = forward1[j];
        if(read1.count > count1){
          count1 = read1.count;
          side1 = 0;
          maxInd1 = j;
        }
      }
      for(int j = 0; j < reverse1.length; j++){
        MappedRead read1 = reverse1[j];
        if(read1.count > count1){
          count1 = read1.count;
          side1 = 1;
          maxInd1 = j;
        }
      }
      for(int k = 0; k < reverse2.length; k++){
        MappedRead read2 = reverse2[k];
        if(read2.count > count2){
          count2 = read2.count;
          side2 = 1;
          maxInd2 = k;
        }
      }
      for(int k = 0; k < forward2.length; k++){
        MappedRead read2 = forward2[k];
        if(read2.count > count2){
          count2 = read2.count;
          side2 = 0;
          maxInd2 = k;
        }
      }
      for(int j = 0; j < forward1.length; j++){
        MappedRead read1 = forward1[j];
        for(int k = 0; k < reverse2.length; k++){
          MappedRead read2 = reverse2[k];
          if(read1.start < read2.end && read2.end - read1.start <= insertLen){
            if(read1.count + read2.count > maxCount){
              maxCount = read1.count + read2.count;
              max1 = j;
              max2 = k;
              side = 0;
            }
          }
        }
      }
      for(int j = 0; j < reverse1.length; j++){
        MappedRead read1 = reverse1[j];
        for(int k = 0; k < forward2.length; k++){
          MappedRead read2 = forward2[k];
          if(read2.start < read1.end && read1.end - read2.start <= insertLen){
            if(read1.count + read2.count > maxCount){
              maxCount = read1.count + read2.count;
              max1 = j;
              max2 = k;
              side = 1;
            }
          }
        }
      }
      switch(side){
        case 0:
          read.r1 = forward1[max1];
          read.r2 = reverse2[max2];
          read.r2.reverse = 1;
          read.r1.q = read.q1;
          read.r2.q = read.q2;
          read.r1.n = 1;
          read.r2.n = 2;
          read.r1.name = read.name;
          read.r2.name = read.name;
          mappedReads.add(read.r1);
          mappedReads.add(read.r2);
          concordant.add(read);
          totalScore += read.r1.count + read.r2.count;
          break;
        case 1:
          read.r1 = reverse1[max1];
          read.r2 = forward2[max2];
          read.r1.reverse = 1;
          read.r1.q = read.q1;
          read.r2.q = read.q2;
          read.r1.n = 1;
          read.r2.n = 2;
          read.r1.name = read.name;
          read.r2.name = read.name;
          mappedReads.add(read.r1);
          mappedReads.add(read.r2);
          concordant.add(read);
          totalScore += read.r1.count + read.r2.count;
          break;
        default:
          if(side1 == -1){
            if(side2 == -1){
              unmapped.add(read);
            }else{
              rightMateMapped.add(read);
              if(side2 == 0){
                read.r2 = forward2[maxInd2];
                mappedForward ++;
              }else{
                read.r2 = reverse2[maxInd2];
                read.r2.reverse = 1;
                mappedReverse ++;
              }
              mappedReads.add(read.r2);
              read.r2.q = read.q2;
              read.r2.n = 2;
              read.r2.name = read.name;
              totalScore += read.r2.count;
            }
          }else{
            if(side2 == -1){
              leftMateMapped.add(read);
              if(side1 == 0){
                read.r1 = forward1[maxInd1];
                mappedForward ++;
              }else{
                read.r1 = reverse1[maxInd1];
                read.r1.reverse = 1;
                mappedReverse ++;
              }
              mappedReads.add(read.r1);
              read.r1.q = read.q1;
              read.r1.n = 1;
              read.r1.name = read.name;
              totalScore += read.r1.count;
            }else{
              discordant.add(read);
              if(side1 == 0){
                read.r1 = forward1[maxInd1];
                mappedForward ++;
              }else{
                read.r1 = reverse1[maxInd1];
                read.r1.reverse = 1;
                mappedReverse ++;
              }
              if(side2 == 0){
                read.r2 = forward2[maxInd2];
                mappedForward ++;
              }else{
                read.r2 = reverse2[maxInd2];
                read.r2.reverse = 1;
                mappedReverse ++;
              }
              mappedReads.add(read.r1);
              mappedReads.add(read.r2);
              read.r1.q = read.q1;
              read.r2.q = read.q2;
              read.r1.n = 1;
              read.r2.n = 2;
              read.r1.name = read.name;
              read.r2.name = read.name;
              totalScore += read.r1.count + read.r2.count;
            }
          }
      }
    }

    @Override
    public void run(){
      mappedReads = new ArrayList<>();
      concordant = new ArrayList<>();
      discordant = new ArrayList<>();
      leftMateMapped = new ArrayList<>();
      rightMateMapped = new ArrayList<>();
      unmapped = new ArrayList<>();
      for(int i = from; i < to; i++){
        PairedRead read = reads.get(i);
        if(!Arrays.equals(read.seq1, NULL_SEQ)){
          if(!Arrays.equals(read.seq2, NULL_SEQ)){
            bothExist++;
            readsTotal += 2;
          }else{
            readsTotal ++;
          }
        }else{
          if(!Arrays.equals(read.seq2, NULL_SEQ)){
            readsTotal ++;
          }
        }
        mapRead(read);
      }
    }

  }

  MappedData mapAllReads(ArrayList<PairedRead> reads, ReferenceAlignment refAlignment) throws InterruptedException{
    long time = System.currentTimeMillis();
    MappedData res = new MappedData();
    int totalPairs = reads.size();
    int totalMapped = 0;
    int mappedForward = 0;
    int mappedReverse = 0;
    int bothExist = 0;
    int bothMapped = 0;
    long totalScore = 0;
    int totalConcordant = 0;
    int readsTotal = 0;
    ArrayList<ReadMapping> tasks = new ArrayList<>();
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    for(int i = 0; i < reads.size(); i += batchSize){
      ReadMapping task = new ReadMapping(reads, refAlignment, i,
          Math.min(i + batchSize, reads.size()));
      executor.execute(task);
      tasks.add(task);
    }
    executor.shutdown();
    executor.awaitTermination(10, TimeUnit.HOURS);
    int mappedReadNum = 0;
    for(ReadMapping task: tasks){
      mappedReadNum += task.mappedReads.size();
    }
    int[] pairs = new int[mappedReadNum];
    byte[] isConcordant = new byte[mappedReadNum];
    int c = 0;
    for(ReadMapping task: tasks){
      for(PairedRead read: task.concordant){
        res.mappedReads.add(read.r1);
        pairs[c] = c + 1;
        isConcordant[c] = 1;
        res.mappedReads.add(read.r2);
        pairs[c + 1] = c;
        isConcordant[c + 1] = 1;
        c += 2;
      }
      for(PairedRead read: task.leftMateMapped){
        res.mappedReads.add(read.r1);
        pairs[c] = -1;
        c++;
      }
      for(PairedRead read: task.rightMateMapped){
        res.mappedReads.add(read.r2);
        pairs[c] = -1;
        c++;
      }
      for(PairedRead read: task.discordant){
        res.mappedReads.add(read.r1);
        pairs[c] = c + 1;
        res.mappedReads.add(read.r2);
        pairs[c + 1] = c;
        c += 2;
      }
      res.concordant.addAll(task.concordant);
      res.discordant.addAll(task.discordant);
      res.leftMateMapped.addAll(task.leftMateMapped);
      res.rightMateMapped.addAll(task.rightMateMapped);
      res.unmapped.addAll(task.unmapped);
      totalMapped += task.mappedReads.size();
      mappedForward += task.mappedForward + task.concordant.size();
      mappedReverse += task.mappedReverse + task.concordant.size();
      totalScore += task.totalScore;
      totalConcordant += task.concordant.size();
      bothExist += task.bothExist;
      bothMapped += task.concordant.size() + task.discordant.size();
      readsTotal += task.readsTotal;
    }
    res.pairs = pairs;
    res.isConcordant = isConcordant;
    logger.println("Mapping to MSA stats:");
    logger.printf("1) Total read pairs: %d\n", totalPairs);
    logger.printf("2) Total pairs with both reads exists: %d\n", bothExist);
    logger.printf("3) Total reads: %d\n", readsTotal);
    logger.printf("4) Total reads mapped from (3): %1.2f, forward: %1.2f, reverse: %1.2f\n",
        (float)totalMapped/readsTotal, (float)mappedForward/totalMapped,
        (float)mappedReverse/totalMapped);
    logger.printf("5) Total pairs with both reads mapped from (2): %d, %1.2f\n",
        bothMapped, (float)bothMapped/bothExist);
    logger.printf("6) Concordant pairs from (5): %1.2f\n", (float)totalConcordant/bothMapped);
    logger.printf("7) Total score: %d, average score: %1.2f\n", totalScore,
        (float) totalScore/totalMapped);
    logger.printf("Total time, s: %d\n", (System.currentTimeMillis() - time)/1000);
    return res;
  }

}
