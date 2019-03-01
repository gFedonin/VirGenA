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
 * Created by Gennady
 */
public class RandomModel extends Constants{

  private KMerCounter counter;
  private Reference genome;
  private int threadNum;
  private int scoreK;
  private int modelK;
  private float coef;
  int[][] counts;

  RandomModel(Document config){
    Element root = config.getRootElement();
    Element elem = root.getChild("Mapper").getChild("RandomModelParameters");
    String pathToReference = root.getChildText("Reference");
    scoreK = Integer.parseInt(root.getChild("Mapper").getChildText("K"));
    modelK = Integer.parseInt(elem.getChildText("Order"));
    coef = Float.parseFloat(elem.getChildText("IndelToleranceThreshold"));
    threadNum = Integer.parseInt(root.getChildText("ThreadNumber"));
    if(threadNum == -1){
      threadNum = Runtime.getRuntime().availableProcessors();
    }
    genome = new Reference(pathToReference, scoreK, threadNum);
    counter = KMerCounter.getInstance(scoreK, coef);
  }

  RandomModel(String pathToGenome, int scoreK, int modelK, float coef, int threadNum){
    genome = new Reference(pathToGenome, scoreK, threadNum);
    counter = KMerCounter.getInstance(scoreK, coef);
    this.threadNum = threadNum;
    this.scoreK = scoreK;
    this.modelK = modelK;
    this.coef = coef;
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
    String[] genomes;
    if(genome.isFragmented){
      genomes = new String[genome.contigEnds.length];
      for(int i = 0, j = 0; i < genome.contigEnds.length; i++){
        genomes[i] = genome.seq.substring(j, genome.contigEnds[i]);
        j = genome.contigEnds[i];
      }
    }else{
      genomes = new String[1];
      genomes[0] = new String(genome.seqB);
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
    String[] genomes;
    if(genome.isFragmented){
      genomes = new String[genome.contigEnds.length];
      for(int i = 0, j = 0; i < genome.contigEnds.length; i++){
        genomes[i] = genome.seq.substring(j, genome.contigEnds[i]);
        j = genome.contigEnds[i];
      }
    }else{
      genomes = new String[1];
      genomes[0] = new String(genome.seqB);
    }
    MarkovModel model = new MarkovModel(genomes, modelK);
    for(int len = maxReadLen, i = size - 1; i >= 0; len -= step, i--){
//      System.out.println("RandomModel len = " + len);
      genRandomScores(model, counts[i], len);
    }
  }

  private class RandomCountsComputing implements Runnable{

    private int len;
    private int[] counts;
    private int from;
    private int to;
    private MarkovModel model;

    RandomCountsComputing(int len, int[] counts, int from, int to, MarkovModel model){
      this.len = len;
      this.counts = counts;
      this.from = from;
      this.to = to;
      this.model = model;
    }

    @Override
    public void run(){
      for(int i = from; i < to; i++){
        byte[] randomRead = model.generate(len);
        counts[i] = counter.computeKMerCount(randomRead, genome.index, 0, genome.length);
      }
    }
  }

  private void genRandomScores(MarkovModel model, int[] counts, int len) throws InterruptedException{
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    int size = counts.length/threadNum;
    for(int i = 0; i < threadNum - 1; i++){
      RandomCountsComputing task = new RandomCountsComputing(len, counts, i*size, i*size + size,
          model);
      executor.execute(task);
    }
    RandomCountsComputing task = new RandomCountsComputing(len, counts, (threadNum - 1)*size,
        counts.length, model);
    executor.execute(task);
    executor.shutdown();
    executor.awaitTermination(10, TimeUnit.HOURS);
    Arrays.sort(counts);
  }

}
