import org.jdom2.Document;
import org.jdom2.Element;

import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Created by Геннадий on 03.03.2015.
 */
public class RandomModelMSA{

  private KMerCounterMSA counter;
  private int threadNum;
  private int scoreK;
  private int modelK;
  private float coef;
  int[][] counts;
  private ReferenceAlignment referenceAlignment;

  RandomModelMSA(Document config){
    Element root = config.getRootElement();
    threadNum = Integer.parseInt(root.getChildText("ThreadNumber"));
    if(threadNum == -1){
      threadNum = Runtime.getRuntime().availableProcessors();
    }
    Element elem = root.getChild("ReferenceSelector");
    String pathToMSA = elem.getChildText("ReferenceMSA");
    scoreK = Integer.parseInt(elem.getChild("MapperToMSA").getChildText("K"));
    elem = elem.getChild("MapperToMSA").getChild("RandomModelParameters");
    modelK = Integer.parseInt(elem.getChildText("Order"));
    coef = Float.parseFloat(elem.getChildText("IndelToleranceThreshold"));
    counter = KMerCounterMSA.getInstance(scoreK, coef);
    referenceAlignment = ReferenceAlignment.getInstance(scoreK, pathToMSA, threadNum);
  }

  RandomModelMSA(String pathToMSA, int scoreK, int modelK, float coef, int threadNum){
    this.threadNum = threadNum;
    this.scoreK = scoreK;
    this.modelK = modelK;
    this.coef = coef;
    counter = KMerCounterMSA.getInstance(scoreK, coef);
    referenceAlignment = ReferenceAlignment.getInstance(scoreK, pathToMSA, threadNum);
  }

  void printModel(int readNum, int minReadLen,
                          int maxReadLen, int step, String outPath)
      throws IOException, InterruptedException{
    DataOutputStream outputStream = new DataOutputStream(new FileOutputStream(outPath));
    outputStream.writeInt(scoreK);
    outputStream.writeFloat(coef);
    outputStream.writeInt(minReadLen);
    outputStream.writeInt(maxReadLen);
    outputStream.writeInt(step);
    outputStream.writeInt(readNum);
    int[] counts = new int[readNum];
    String[] genomes = new String[referenceAlignment.refAlns.size()];
    for(int i = 0; i < genomes.length; i++){
      genomes[i] = referenceAlignment.refAlns.get(i).seq;
    }
    MarkovModel model = new MarkovModel(genomes, modelK);
    for(int len = minReadLen; len <= maxReadLen; len += step){
      System.out.println("RandomModel len = " + len);
      genRandomScores(model, counts, len);
      for(int i = 0; i < readNum; i++){
        outputStream.writeInt(counts[i]);
      }
    }
    outputStream.close();
  }

  void genModel(int readNum, int minReadLen,
                  int maxReadLen, int step)
      throws InterruptedException{
    int size = (maxReadLen - minReadLen)/step + 1;
    counts = new int[size][readNum];
    String[] genomes = new String[referenceAlignment.refAlns.size()];
    for(int i = 0; i < genomes.length; i++){
      genomes[i] = referenceAlignment.refAlns.get(i).seq;
    }
    MarkovModel model = new MarkovModel(genomes, modelK);
    for(int len = maxReadLen, i = size - 1; i >= 0; len -= step, i--){
//      System.out.println("RandomModel len = " + len);
      genRandomScores(model, counts[i], len);
    }
  }

  private class RandomCountsComputing implements Runnable{

    public int len;
    public int[] counts;
    public int from;
    public int to;
    public MarkovModel model;
    private ReferenceAlignment referenceAlignment;

    public RandomCountsComputing(int len, int[] counts, int from, int to, MarkovModel model,
                                 ReferenceAlignment referenceAlignment){
      this.len = len;
      this.counts = counts;
      this.from = from;
      this.to = to;
      this.model = model;
      this.referenceAlignment = referenceAlignment;
    }

    @Override
    public void run(){
      for(int i = from; i < to; i++){
        byte[] randomRead = model.generate(len);
        counts[i] = counter.computeKMerCount(randomRead, referenceAlignment.refAlnIndex,
            0, referenceAlignment.length);
      }
    }
  }

  private void genRandomScores(MarkovModel model, int[] counts, int len) throws InterruptedException{
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    int size = counts.length/threadNum;
    for(int i = 0; i < threadNum - 1; i++){
      RandomCountsComputing task = new RandomCountsComputing(len, counts, i*size, i*size + size, model, referenceAlignment);
      executor.execute(task);
    }
    RandomCountsComputing task = new RandomCountsComputing(len, counts, (threadNum - 1)*size, counts.length, model, referenceAlignment);
    executor.execute(task);
    executor.shutdown();
    executor.awaitTermination(10, TimeUnit.HOURS);
    Arrays.sort(counts);
  }

}
