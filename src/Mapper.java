import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentGroup;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;
import org.jdom2.Document;
import org.jdom2.input.SAXBuilder;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Created by Gennady
 */
public class Mapper extends Constants{

  private int insertLen;
  KMerCounter counter;
  Aligner aligner;
  Logger logger;
  private int threadNum;
  private int batchSize;

  Mapper(Document document){
    aligner = new SmithWatermanGotoh(document);
    counter = KMerCounter.getInstance(document);
    logger = Logger.getInstance(document);
    insertLen = Integer.parseInt(document.getRootElement().getChild("Data").getChildText("InsertionLength"));
    threadNum = Integer.parseInt(document.getRootElement().getChildText("ThreadNumber"));
    if(threadNum == -1){
      threadNum = Runtime.getRuntime().availableProcessors();
    }
    batchSize = Integer.parseInt(document.getRootElement().getChildText("BatchSize"));
  }

  Mapper(Document document, Namespace parsedArgs){
    aligner = new SmithWatermanGotoh(document);
    counter = KMerCounter.getInstance(document);
    logger = Logger.getInstance(document);
    insertLen = Integer.parseInt(document.getRootElement().getChild("Data").getChildText("InsertionLength"));
    String n = parsedArgs.getString("thread_num");
    if(n == null){
      threadNum = Integer.parseInt(document.getRootElement().getChildText("ThreadNumber"));
    }else{
      threadNum = Integer.parseInt(n);
    }
    if(threadNum == -1){
      threadNum = Runtime.getRuntime().availableProcessors();
    }
    batchSize = Integer.parseInt(document.getRootElement().getChildText("BatchSize"));
  }

  private MappedRead[] mapReadFast(byte[] read, Reference genome){
    MappedRead[] bestPositions = counter.getNBestRegions(read, genome.index, genome.contigEnds);
    if(bestPositions == null){
      return new MappedRead[0];
    }
    for(MappedRead mappedRead : bestPositions){
      mappedRead.seq = read;
    }
    return bestPositions;
  }

  private class ReadsMapping implements Runnable{

    private ArrayList<PairedRead> reads;
    private int from;
    private int to;
    ArrayList<MappedRead> mappedReads;
    ArrayList<PairedRead> concordant;
    ArrayList<PairedRead> discordant;
    ArrayList<PairedRead> leftMateMapped;
    ArrayList<PairedRead> rightMateMapped;
    ArrayList<PairedRead> unmapped;
    int readsTotal = 0;
    int bothExist = 0;
    int mappedForward = 0;
    int mappedReverse = 0;
    long totalScore = 0;
    private Reference genome;

    ReadsMapping(ArrayList<PairedRead> reads, int from, int to, Reference genome){
      this.reads = reads;
      this.from = from;
      this.to = to;
      this.genome = genome;
    }

    private void mapRead(PairedRead read){
//      if(read.name.equals("700150F:34:HVJVGBCXY:1:1103:17680:52747")){
//        int a = 0;
//      }
      MappedRead[] forward1 = mapReadFast(read.seq1, genome);
      MappedRead[] reverse1 = mapReadFast(getComplement(read.seq1), genome);
      MappedRead[] forward2 = mapReadFast(read.seq2, genome);
      MappedRead[] reverse2 = mapReadFast(getComplement(read.seq2), genome);
      int maxCount = 0;
      int max1 = -1;
      int max2 = -1;
      byte side = -1;
      int maxInd1 = -1;
      int count1 = 0;
      int side1 = -1;
      int[] maxIndFrag1 = null;
      int[] countFrag1 = null;
      int[] sideFrag1 = null;
      int maxInd2 = -1;
      int count2 = 0;
      int side2 = -1;
      int[] maxIndFrag2 = null;
      int[] countFrag2 = null;
      int[] sideFrag2 = null;
      if(genome.isFragmented){
        maxIndFrag1 = new int[genome.fragmentEnds.length];
        countFrag1 = new int[genome.fragmentEnds.length];
        sideFrag1 = new int[genome.fragmentEnds.length];
        maxIndFrag2 = new int[genome.fragmentEnds.length];
        countFrag2 = new int[genome.fragmentEnds.length];
        sideFrag2 = new int[genome.fragmentEnds.length];
      }
      for(int j = 0; j < forward1.length; j++){
        MappedRead read1 = forward1[j];
        if(genome.isFragmented){
          read1.chooseFragment(genome.fragmentEnds);
          int fragmentID1 = read1.fragmentID;
          if(read1.count > countFrag1[fragmentID1]){
            countFrag1[fragmentID1] = read1.count;
            sideFrag1[fragmentID1] = 0;
            maxIndFrag1[fragmentID1] = j;
          }
        }
        if(read1.count > count1){
          count1 = read1.count;
          side1 = 0;
          maxInd1 = j;
        }
      }
      for(int j = 0; j < reverse1.length; j++){
        MappedRead read1 = reverse1[j];
        if(genome.isFragmented){
          read1.chooseFragment(genome.fragmentEnds);
          int fragmentID1 = read1.fragmentID;
          if(read1.count > countFrag1[fragmentID1]){
            countFrag1[fragmentID1] = read1.count;
            sideFrag1[fragmentID1] = 1;
            maxIndFrag1[fragmentID1] = j;
          }
        }
        if(read1.count > count1){
          count1 = read1.count;
          side1 = 1;
          maxInd1 = j;
        }
      }
      for(int k = 0; k < reverse2.length; k++){
        MappedRead read2 = reverse2[k];
        if(genome.isFragmented){
          read2.chooseFragment(genome.fragmentEnds);
          int fragmentID2 = read2.fragmentID;
          if(read2.count > countFrag2[fragmentID2]){
            countFrag2[fragmentID2] = read2.count;
            sideFrag2[fragmentID2] = 1;
            maxIndFrag2[fragmentID2] = k;
          }
        }
        if(read2.count > count2){
          count2 = read2.count;
          side2 = 1;
          maxInd2 = k;
        }
      }
      for(int k = 0; k < forward2.length; k++){
        MappedRead read2 = forward2[k];
        if(genome.isFragmented){
          read2.chooseFragment(genome.fragmentEnds);
          int fragmentID2 = read2.fragmentID;
          if(read2.count > countFrag2[fragmentID2]){
            countFrag2[fragmentID2] = read2.count;
            sideFrag2[fragmentID2] = 0;
            maxIndFrag2[fragmentID2] = k;
          }
        }
        if(read2.count > count2){
          count2 = read2.count;
          side2 = 0;
          maxInd2 = k;
        }
      }

      for(int j = 0; j < forward1.length; j++){
        MappedRead read1 = forward1[j];
        int fragmentID1 = read1.fragmentID;
        for(int k = 0; k < reverse2.length; k++){
          MappedRead read2 = reverse2[k];
          int fragmentID2 = read2.fragmentID;
          if(read1.start < read2.end && read2.end - read1.start <= insertLen){
            if(fragmentID1 == fragmentID2 && read1.count + read2.count > maxCount){
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
        int fragmentID1 = read1.fragmentID;
        for(int k = 0; k < forward2.length; k++){
          MappedRead read2 = forward2[k];
          int fragmentID2 = read2.fragmentID;
          if(read2.start < read1.end && read1.end - read2.start <= insertLen){
            if(fragmentID1 == fragmentID2 && read1.count + read2.count > maxCount){
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
          byte[] genomeSubseq = Arrays.copyOfRange(genome.seqB, read.r1.start, read.r1.end);
          read.r1.aln = aligner.align(read.r1.seq, genomeSubseq);
          genomeSubseq = Arrays.copyOfRange(genome.seqB, read.r2.start, read.r2.end);
          read.r2.aln = aligner.align(read.r2.seq, genomeSubseq);
          mappedReads.add(read.r1);
          mappedReads.add(read.r2);
          concordant.add(read);
          totalScore += read.r1.aln.score + read.r2.aln.score;
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
          genomeSubseq = Arrays.copyOfRange(genome.seqB, read.r1.start, read.r1.end);
          read.r1.aln = aligner.align(read.r1.seq, genomeSubseq);
          genomeSubseq = Arrays.copyOfRange(genome.seqB, read.r2.start, read.r2.end);
          read.r2.aln = aligner.align(read.r2.seq, genomeSubseq);
          mappedReads.add(read.r1);
          mappedReads.add(read.r2);
          concordant.add(read);
          totalScore += read.r1.aln.score + read.r2.aln.score;
          break;
        default:
          if(side1 == -1){
            if(side2 == -1){
              unmapped.add(read);
              read.r1 = null;
              read.r2 = null;
            }else{
              rightMateMapped.add(read);
              read.r1 = null;
              if(side2 == 0){
                read.r2 = forward2[maxInd2];
                mappedForward ++;
              }else{
                read.r2 = reverse2[maxInd2];
                read.r2.reverse = 1;
                mappedReverse ++;
              }
              genomeSubseq = Arrays.copyOfRange(genome.seqB, read.r2.start, read.r2.end);
              read.r2.aln = aligner.align(read.r2.seq, genomeSubseq);
              mappedReads.add(read.r2);
              read.r2.q = read.q2;
              read.r2.n = 2;
              read.r2.name = read.name;
              totalScore += read.r2.aln.score;
            }
          }else{
            if(side2 == -1){
              leftMateMapped.add(read);
              read.r2 = null;
              if(side1 == 0){
                read.r1 = forward1[maxInd1];
                mappedForward ++;
              }else{
                read.r1 = reverse1[maxInd1];
                read.r1.reverse = 1;
                mappedReverse ++;
              }
              genomeSubseq = Arrays.copyOfRange(genome.seqB, read.r1.start, read.r1.end);
              read.r1.aln = aligner.align(read.r1.seq, genomeSubseq);
              mappedReads.add(read.r1);
              read.r1.q = read.q1;
              read.r1.n = 1;
              read.r1.name = read.name;
              totalScore += read.r1.aln.score;
            }else{
              if(genome.isFragmented){
                int maxSum = 0;
                int maxIndex = -1;
                for(int i = 0; i < genome.fragmentEnds.length; i++){
                  if(countFrag1[i] > 0 && countFrag2[i] > 0 && countFrag1[i] + countFrag2[i] > maxSum){
                    maxSum = countFrag1[i] + countFrag2[i];
                    maxIndex = i;
                  }
                }
                if(maxIndex != -1){
                  discordant.add(read);
                  if(sideFrag1[maxIndex] == 0){
                    read.r1 = forward1[maxIndFrag1[maxIndex]];
                    mappedForward ++;
                  }else{
                    read.r1 = reverse1[maxIndFrag1[maxIndex]];
                    read.r1.reverse = 1;
                    mappedReverse ++;
                  }
                  if(sideFrag2[maxIndex] == 0){
                    read.r2 = forward2[maxIndFrag2[maxIndex]];
                    mappedForward ++;
                  }else{
                    read.r2 = reverse2[maxIndFrag2[maxIndex]];
                    read.r2.reverse = 1;
                    mappedReverse ++;
                  }
                  genomeSubseq = Arrays.copyOfRange(genome.seqB, read.r1.start, read.r1.end);
                  read.r1.aln = aligner.align(read.r1.seq, genomeSubseq);
                  genomeSubseq = Arrays.copyOfRange(genome.seqB, read.r2.start, read.r2.end);
                  read.r2.aln = aligner.align(read.r2.seq, genomeSubseq);
                  mappedReads.add(read.r1);
                  mappedReads.add(read.r2);
                  read.r1.q = read.q1;
                  read.r2.q = read.q2;
                  read.r1.n = 1;
                  read.r2.n = 2;
                  read.r1.name = read.name;
                  read.r2.name = read.name;
                  totalScore += read.r1.aln.score + read.r2.aln.score;
                }else{
                  unmapped.add(read);
                  read.r1 = null;
                  read.r2 = null;
                }
              }else{
                discordant.add(read);
                if(side1 == 0){
                  read.r1 = forward1[maxInd1];
                  mappedForward++;
                }else{
                  read.r1 = reverse1[maxInd1];
                  read.r1.reverse = 1;
                  mappedReverse++;
                }
                if(side2 == 0){
                  read.r2 = forward2[maxInd2];
                  mappedForward++;
                }else{
                  read.r2 = reverse2[maxInd2];
                  read.r2.reverse = 1;
                  mappedReverse++;
                }
                genomeSubseq = Arrays.copyOfRange(genome.seqB, read.r1.start, read.r1.end);
                read.r1.aln = aligner.align(read.r1.seq, genomeSubseq);
                genomeSubseq = Arrays.copyOfRange(genome.seqB, read.r2.start, read.r2.end);
                read.r2.aln = aligner.align(read.r2.seq, genomeSubseq);
                mappedReads.add(read.r1);
                mappedReads.add(read.r2);
                read.r1.q = read.q1;
                read.r2.q = read.q2;
                read.r1.n = 1;
                read.r2.n = 2;
                read.r1.name = read.name;
                read.r2.name = read.name;
                totalScore += read.r1.aln.score + read.r2.aln.score;
              }

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

  MappedData mapReads(ArrayList<PairedRead> reads, Reference genome) throws InterruptedException{
    MappedData res = new MappedData();
    if(genome.length == 0){
      res.unmapped = reads;
      logger.println("Mapping is impossible: reference is an empty string!");
      return res;
    }
    int totalPairs = reads.size();
    int totalMapped = 0;
    int mappedForward = 0;
    int mappedReverse = 0;
    int bothExist = 0;
    int bothMapped = 0;
    long totalScore = 0;
    int totalConcordant = 0;
    int readsTotal = 0;
    ArrayList<ReadsMapping> tasks = new ArrayList<>();
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    for(int i = 0; i < reads.size(); i += batchSize){
      ReadsMapping task = new ReadsMapping(reads, i,
          Math.min(i + batchSize, reads.size()), genome);
      executor.execute(task);
      tasks.add(task);
    }
    executor.shutdown();
    executor.awaitTermination(10, TimeUnit.HOURS);
    for(ReadsMapping task: tasks){
      res.mappedReads.addAll(task.mappedReads);
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
    logger.println("Mapping stats:");
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
        (float) totalScore / totalMapped);
    return res;
  }

  void mapReads(MappedData mappedData, Reference genome) throws InterruptedException{
    logger.printf("Remapping %d read pairs\n", mappedData.needToRemap.size());
    int totalMapped = mappedData.mappedReads.size();
    int mappedForward = 0;
    int mappedReverse = 0;
    int bothExist = mappedData.concordant.size();
    int bothMapped = mappedData.concordant.size();
    long totalScore = 0;
    int totalConcordant = mappedData.concordant.size();
    int readsTotal = mappedData.mappedReads.size();
    for(MappedRead mappedRead: mappedData.mappedReads){
      if(mappedRead.reverse == 0){
        mappedForward ++;
      }else{
        mappedReverse ++;
      }
      totalScore += mappedRead.aln.score;
    }
    int totalPairs = mappedData.needToRemap.size() + mappedData.concordant.size() +
        mappedData.leftMateMapped.size() + mappedData.rightMateMapped.size();
    ArrayList<ReadsMapping> tasks = new ArrayList<>();
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    for(int i = 0; i < mappedData.needToRemap.size(); i += batchSize){
      ReadsMapping task = new ReadsMapping(mappedData.needToRemap, i,
          Math.min(i + batchSize, mappedData.needToRemap.size()), genome);
      executor.execute(task);
      tasks.add(task);
    }
    executor.shutdown();
    executor.awaitTermination(10, TimeUnit.HOURS);
    for(ReadsMapping task: tasks){
      mappedData.mappedReads.addAll(task.mappedReads);
      mappedData.concordant.addAll(task.concordant);
      mappedData.discordant.addAll(task.discordant);
      mappedData.leftMateMapped.addAll(task.leftMateMapped);
      mappedData.rightMateMapped.addAll(task.rightMateMapped);
      mappedData.unmapped.addAll(task.unmapped);
      totalMapped += task.mappedReads.size();
      mappedForward += task.mappedForward + task.concordant.size();
      mappedReverse += task.mappedReverse + task.concordant.size();
      totalScore += task.totalScore;
      totalConcordant += task.concordant.size();
      bothExist += task.bothExist;
      bothMapped += task.concordant.size() + task.discordant.size();
      readsTotal += task.readsTotal;
    }
    mappedData.needToRemap.clear();
    logger.println("Mapping stats:");
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
        (float) totalScore / totalMapped);
  }

  static void addParameters(ArgumentParser parser){
    parser.description("This tool will map paired reads " +
        "to given reference sequence. The resulting sorted and indexed BAM file 'mapped_reads.bam' will be created in " +
        "the given output folder along with 'mapped_reads.bai' and log file.");
    parser.addArgument("-c").dest("config_file").help("Path to config file").required(true);
    ArgumentGroup fileGroup = parser.addArgumentGroup("Input files").description("Paths to fastq files with reads to map. " +
        "May be gzipped. If present, the same paths in config file are ignored.");
    fileGroup.addArgument("-p1").dest("pair1").help("Path to the first paired read file or a list of comma separated files");
    fileGroup.addArgument("-p2").dest("pair2").help("Path to the second paired read file or a list of comma separated files");
    parser.addArgument("-r").dest("reference").help("Path to the reference FASTA file");
    parser.addArgument("-o").dest("output_path").help("The output folder. Intermediate folders will be " +
        "created if not exists. Of present, the output path in config file will be ignored. If not - 'mapped_reads.bam' will be " +
        "created in the output dir from config file. Default output dir is current dir.");
    parser.addArgument("-n").dest("thread_num").help("Number of threads to use, set to -1 to use all available processors. " +
        "If not specified - the value from config file will be used.");
  }

  static void run(Namespace parsedArgs){
    try{
//      Namespace parsedArgs = parser.parseArgs(args);
      SAXBuilder jdomBuilder = new SAXBuilder();
      Document jdomDocument = jdomBuilder.build(parsedArgs.getString("config_file"));
      long time = System.currentTimeMillis();
      DataReader dataReader = DataReader.getInstance(jdomDocument, parsedArgs);
      Mapper mapper = new Mapper(jdomDocument, parsedArgs);
      Reference genome = new Reference(jdomDocument, parsedArgs);
      MappedData mappedData = mapper.mapReads(dataReader.pairedReads, genome);
      mapper.logger.printf("Time: %d, s\n", (System.currentTimeMillis() - time) / 1000);
      BamPrinter bamPrinter = new BamPrinter();
      String outPath = parsedArgs.getString("output_dir");
      if(outPath == null){
        outPath = jdomDocument.getRootElement().getChildText("OutPath");
        if(outPath.isEmpty()){
          outPath = System.getProperty("user.dir");
        }
      }
      bamPrinter.printBAM(mappedData, genome, outPath + "mapped_reads.bam");
    }catch(Exception e){
      e.printStackTrace();
    }
  }

  public static void main(String[] args){
    ArgumentParser parser = ArgumentParsers.newFor("Mapper").build();
    addParameters(parser);
    try{
      Namespace parsedArgs = parser.parseArgs(args);
      run(parsedArgs);
    }catch(ArgumentParserException e){
      parser.handleError(e);
    }
  }

}
