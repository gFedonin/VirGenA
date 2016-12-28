import java.io.*;
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
    for(int len = minReadLen; len <= maxReadLen; len += step){
      System.out.println("RandomModel len = " + len);
      genRandomScores(counts, len);
      for(int i = 0; i < readNum; i++){
        outputStream.writeInt(counts[i]);
      }
    }
    outputStream.close();
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

  private void genRandomScores(int[] counts, int len) throws InterruptedException{
    RandomCountsComputing[] tasks = new RandomCountsComputing[threadNum];
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    int size = counts.length/threadNum;
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
    for(int i = 0; i < threadNum - 1; i++){
      RandomCountsComputing task = new RandomCountsComputing(len, counts, i*size, i*size + size,
          model);
      executor.execute(task);
      tasks[i] = task;
    }
    RandomCountsComputing task = new RandomCountsComputing(len, counts, (threadNum - 1)*size,
        counts.length, model);
    executor.execute(task);
    tasks[threadNum - 1] = task;
    executor.shutdown();
    executor.awaitTermination(10, TimeUnit.HOURS);
    Arrays.sort(counts);
  }

}
