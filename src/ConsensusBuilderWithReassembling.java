import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.jdom2.Document;
import org.jdom2.Element;
import org.jdom2.input.SAXBuilder;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.StringTokenizer;

/**
 * Created by Геннадий on 06.08.2014.
 */
public class ConsensusBuilderWithReassembling extends ConsensusBuilderSimple{

  private ArrayList<ProblemRead> problemReads;
  protected ArrayList<PairedRead> reads;
  private ProblemReadMapper problemReadMapper;
  Interval[] problemIntervals;
  private byte[] simpleConsensus;
  private float identityThreshold;
  private int minLen;
  private int minFinalReads;
  int K;
  private int minIntersectLen;
  Logger logger;
  private String outPath;
  String finalConsensus;
  int[] contigEnds;
  int[] fragmentEnds;
  String[] contigs;
  ArrayList<byte[]>[] contigsFragmented;
  boolean genomeIsFragmented;
  String[] fragmentNames;
  int[] genomeToAssembly;
  String genomeName;
  private boolean debug;

  ConsensusBuilderWithReassembling(Document document){
    super(document);
    Element cbElement = document.getRootElement().getChild("ConsensusBuilder");
    Element raElement = cbElement.getChild("Reassembler");
    debug = Boolean.parseBoolean(cbElement.getChildText("Debug"));
    identityThreshold = Float.parseFloat(cbElement.getChildText("IdentityThreshold"));
    minLen = Integer.parseInt(raElement.getChildText("MinReadLength"));
    minFinalReads = Integer.parseInt(cbElement.getChildText("MinTerminationReadsNumber"));
    K = KMerCounter.getInstance(document).K;
    minIntersectLen = Integer.parseInt(cbElement.getChildText("MinIntersectionLength"));
    logger = Logger.getInstance(document);
    outPath = document.getRootElement().getChildText("OutPath");
    problemReadMapper = new ProblemReadMapperPairedReadTermination(document);
  }

  private ProblemRead assignToGap(MappedRead read){
    if(read == null || read.seq.length < minLen){
      return null;
    }
    Interval temp = new Interval();
    temp.start = (short)(read.aln.start2 + read.start);
    temp.end = (short)(read.aln.end2 + read.start);
    short i = (short)Arrays.binarySearch(problemIntervals, temp);
    if(i < 0){
      i = (short)(-i - 1);
    }
    int readLen = read.seq.length;
    short cutoff = (short)mapper.counter.cutoffs[readLen];
    if(i == problemIntervals.length){
      Interval prev = problemIntervals[i - 1];
      if(temp.start - read.aln.start1 < prev.end){
        ProblemRead pRead = new ProblemRead(read.name, read.n, read.seq, (short)(i - 1), read.reverse);
        pRead.count = -1;
        pRead.cutoff = cutoff;
        problemReads.add(pRead);
        return pRead;
      }
      return null;
    }
    if(i == 0){
      Interval next = problemIntervals[i];
      if(temp.end + read.seq.length - read.aln.end1 > next.start){
        ProblemRead pRead = new ProblemRead(read.name, read.n, read.seq, i, read.reverse);
        pRead.count = -1;
        pRead.cutoff = cutoff;
        problemReads.add(pRead);
        return pRead;
      }
      return null;
    }
    Interval prev = problemIntervals[i - 1];
    Interval next = problemIntervals[i];
    int intersectsPrev = prev.end - temp.start + read.aln.start1;
    int intersectsNext = temp.end + read.seq.length - read.aln.end1 - next.start;
    if(intersectsNext > 0){
      if(intersectsPrev > intersectsNext){
        ProblemRead pRead = new ProblemRead(read.name, read.n, read.seq, (short)(i - 1), read.reverse);
        pRead.count = -1;
        pRead.cutoff = cutoff;
        problemReads.add(pRead);
        return pRead;
      }
      ProblemRead pRead = new ProblemRead(read.name, read.n, read.seq, i, read.reverse);
      pRead.count = -1;
      pRead.cutoff = cutoff;
      problemReads.add(pRead);
      return  pRead;
    }else if(intersectsPrev > 0){
      ProblemRead pRead = new ProblemRead(read.name, read.n, read.seq, (short)(i - 1), read.reverse);
      pRead.count = -1;
      pRead.cutoff = cutoff;
      problemReads.add(pRead);
      return pRead;
    }
    return null;
  }

  private void selectGapIntersectingReads(ArrayList<PairedRead> reads){
    problemReads = new ArrayList<>();
    for(PairedRead pairedRead: reads){
      ProblemRead pRead1 = assignToGap(pairedRead.r1);
      ProblemRead pRead2 = assignToGap(pairedRead.r2);
      if(pRead1 != null){
        if(pRead2 != null){
          pRead1.mate = pRead2;
        }else{
          pRead1.mate = pairedRead.r2;
        }
      }
      if(pRead2 != null){
        if(pRead1 != null){
          pRead2.mate = pRead1;
        }else{
          pRead2.mate = pairedRead.r1;
        }
      }
    }
  }

  private boolean growLeftEdges() throws InterruptedException{
    boolean someIntervalWasUpdated = false;
    for(Interval gap : problemIntervals){
      if(gap.growLeftEdge(coverageThreshold)){
        someIntervalWasUpdated = true;
        gap.buildIndexLeft(K, threadNum);
        gap.buildIndexConcat(K, threadNum);
      }
    }
    return someIntervalWasUpdated;
  }

  private boolean growRightEdges() throws InterruptedException{
    boolean someIntervalWasUpdated = false;
    for(Interval gap : problemIntervals){
      if(gap.growRightEdge(coverageThreshold)){
        someIntervalWasUpdated = true;
        gap.buildIndexRight(K, threadNum);
        gap.buildIndexConcat(K, threadNum);
      }
    }
    return someIntervalWasUpdated;
  }

  private int estimateAnchorSize(ArrayList<PairedRead> reads){
    int mReadLen = 0;
    for(PairedRead read: reads){
      int len = read.seq1.length;
      if(len > mReadLen){
        mReadLen = len;
      }
      len = read.seq2.length;
      if(len > mReadLen){
        mReadLen = len;
      }
    }
    logger.println("Max read len = " + mReadLen);
    float maxCutoff = mapper.counter.cutoffs[mReadLen];//Math.max(cutoffs[mReadLen], complimentCutoffs[mReadLen]);
    logger.println("K-mer treshold = " + maxCutoff);
    int anchorSize;
    if((1.0f - identityThreshold)*K < 1.0f){
      anchorSize = (int)Math.ceil((maxCutoff + K)/(1.0f - (1.0f - identityThreshold)*K));
    }else{
      anchorSize = 100;
    }
    logger.println("Anchor size = " + anchorSize);
    return anchorSize;
  }

  private void printGaps(){
    int[] gapReadCounts = new int[problemIntervals.length];
    int[] gapReadCountsF = new int[problemIntervals.length];
    int[] gapReadCountsR = new int[problemIntervals.length];
    for(ProblemRead read: problemReads){
      if(read.gapID != -1){
        gapReadCounts[read.gapID] ++;
        if(read.reverse == 0){
          gapReadCountsF[read.gapID] ++;
        }else{
          gapReadCountsR[read.gapID] ++;
        }
      }
    }
    logger.println("Read counts in gaps:");
    for(int i = 0; i < problemIntervals.length; i++){
      Interval gap = problemIntervals[i];
      logger.println(gap.start + " " + gap.end + ": " + gapReadCounts[i] + " F: " +
          (float)gapReadCountsF[i]/gapReadCounts[i] + " R: " +
          (float)gapReadCountsR[i]/gapReadCounts[i]);
    }
  }

  private void mergeEdges(Interval interval, Aligner aligner, float identityThreshold,
                          int minFinalReads, int minIntersectLen){
    if(interval.concat.length == 0){
      interval.seqFinal = "------";
    }
    if(interval.left.length == 0){
      int rightAddedLen = interval.concat.length - anchorSize;
      interval.seqFinal = "------" + new String(interval.concat, 0, rightAddedLen);
      if(debug){
        logger.println("Merging gap: " + interval.start + " " + interval.end);
        logger.println("Left anchor seqB: ");
        logger.println("Left seqB: ");
        logger.println("Right anchor seqB: " + new String(interval.concat,
            rightAddedLen, anchorSize));
        logger.println("Right seqB: " + new String(interval.concat, 0, rightAddedLen));
        logger.println("Seq to insert: " + interval.seqFinal);
      }
      return;
    }
    if(interval.right.length == 0){
      int leftAddedLen = interval.concat.length - anchorSize;
      interval.seqFinal = new String(interval.concat, anchorSize, leftAddedLen) + "------";
      if(debug){
        logger.println("Merging gap: " + interval.start + " " + interval.end);
        logger.println("Left anchor seqB: " + new String(interval.concat, 0, anchorSize));
        logger.println("Left seqB: " + new String(interval.concat, anchorSize, leftAddedLen));
        logger.println("Right anchor seqB: ");
        logger.println("Right seqB: ");
        logger.println("Seq to insert: " + interval.seqFinal);
      }
      return;
    }
    int rightAddedLen = interval.concat.length - interval.junction - anchorSize - 1;
    int leftAddedLen = interval.junction - anchorSize + 1;
    if(debug){
      logger.println("Merging gap: " + interval.start + " " + interval.end);
      logger.println("Found " + interval.finalizingReadsFound + " finalizing reads");
      logger.println("Left anchor seqB: " + new String(interval.concat, 0, anchorSize));
      logger.println("Left seqB: " + new String(interval.concat, anchorSize, leftAddedLen));
      logger.println("Right anchor seqB: " + new String(interval.concat,
          interval.junction + rightAddedLen + 1, anchorSize));
      logger.println("Right seqB: " + new String(interval.concat,
          interval.junction + 1, rightAddedLen));
    }
    if(interval.finalizingReadsFound >= minFinalReads){
      if(debug){
        logger.println("Successfully merged both seqs!");
      }
      interval.seqFinal = new String(interval.concat, anchorSize, leftAddedLen + rightAddedLen);
    }else{
      Alignment aln = aligner.align(interval.left, interval.right);
      int minIntersection = Math.min(aln.end1 - aln.start1, aln.end2 - aln.start2);
      if(minIntersection < minIntersectLen){
        if(debug){
          logger.println("Too small intersection of left and right seqs!");
        }
        interval.seqFinal = new String(interval.concat, anchorSize, leftAddedLen) + "------" +
                new String(interval.concat, interval.junction + 1, rightAddedLen);
      }else if(aln.end1 < anchorSize){
        if(debug){
          logger.println("Left and right seqs intersect inside of left anchor!");
        }
        interval.seqFinal = new String(interval.concat, anchorSize, leftAddedLen) + "------" +
                new String(interval.concat, interval.junction + 1, rightAddedLen);
      }else if(aln.start2 > rightAddedLen){
        if(debug){
          logger.println("Left and right seqs intersect inside of right anchor!");
        }
        interval.seqFinal = new String(interval.concat, anchorSize, leftAddedLen) + "------" +
                new String(interval.concat, interval.junction + 1, rightAddedLen);
      }else if((float)aln.identity/aln.length > identityThreshold){
        if(debug){
          logger.println("Left and right seqs were merged after SW alignment!");
          logger.println("Alignment: left: " + aln.start1 + " - " + aln.end1 + " right: " +
              aln.start2 + " - " + aln.end2);
        }
        byte[] finSeq;
        if(interval.left.length > interval.right.length){
          finSeq  = new byte[aln.end1 + interval.right.length - aln.end2];
          System.arraycopy(interval.left, 0, finSeq, 0, aln.end1);
          System.arraycopy(interval.right, aln.end2, finSeq, aln.end1,
                  interval.right.length - aln.end2);
        }else{
          finSeq  = new byte[aln.start1 + interval.right.length - aln.start2];
          System.arraycopy(interval.left, 0, finSeq, 0, aln.start1);
          System.arraycopy(interval.right, aln.start2, finSeq, aln.start1,
                  interval.right.length - aln.start2);
        }
        if(finSeq.length > 2*anchorSize){
          interval.seqFinal = new String(finSeq, anchorSize, finSeq.length - 2*anchorSize);
        }else{
          interval.seqFinal = new String(interval.concat, anchorSize, leftAddedLen) + "------" +
                  new String(interval.concat, interval.junction + 1, rightAddedLen);
        }
      }else{
        if(debug){
          logger.println("Too small identity to merge with SW!");
        }
        interval.seqFinal = new String(interval.concat, anchorSize, leftAddedLen) + "------" +
                new String(interval.concat, interval.junction + 1, rightAddedLen);
      }
    }
    if(debug){
      logger.println("Seq to insert: " + interval.seqFinal);
    }
  }

//  void indexBam(String bamPath) throws IOException, InterruptedException{
//    Process p = Runtime.getRuntime().exec("samtools index " + bamPath);
//    BufferedReader streamReader = new BufferedReader(new InputStreamReader(p.getInputStream()));
//    String line;
//    while((line = streamReader.readLine()) != null){
//      //System.out.println(line);
//    }
//    streamReader = new BufferedReader(new InputStreamReader(p.getErrorStream()));
//    while((line = streamReader.readLine()) != null){
//      //System.out.println(line);
//    }
//    p.waitFor();
//  }

  @Override
  public String buildConsensus(Reference genome, ArrayList<PairedRead> reads){
    try{
      this.reads = reads;
      genomeIsFragmented = genome.isFragmented;
      fragmentNames = genome.fragmentNames;
      contigEnds = genome.contigEnds;
      fragmentEnds = genome.fragmentEnds;
      genomeName = genome.name;
      int windowSize = (int) Math.ceil(1.0f/(1.0f - identityThreshold));
      //building first order simple assembly
      logger.println("Iteration 1");
      long time = System.currentTimeMillis();
      mappedData = mapper.mapReads(reads, genome);
      if(debug){
        logger.println("Simple consensus:");
        logger.println(new String(getConsensusSimple(false, coverageThreshold, genome.seqB)));
      }
      simpleConsensus = getConsensusSimple(true, coverageThreshold, genome.seqB);
      if(debug){
        printAssembly("consensus1.fasta");
        BamPrinter bamPrinter = new BamPrinter();
        String bamPath = outPath + "mapped_reads1.bam";
        Reference ref1 = new Reference(simpleConsensus, K, genome, false, threadNum);
        bamPrinter.printBAM(mappedData, ref1, bamPath);
      }
      logger.println("1st order consensus with filled gaps:");
      logger.println(new String(simpleConsensus));
      logger.printf("Time iter 1, s: %d\n", (System.currentTimeMillis() - time)/1000);
      logger.println("Iteration 2");
      time = System.currentTimeMillis();
      //change genome to assembly, rebuild indexes and remap all reads to simple assembly
      genome = new Reference(simpleConsensus, K, genome, threadNum);
      mappedData = mapper.mapReads(reads, genome);
      simpleConsensus = getConsensusSimple(false, coverageThreshold, genome.seqB);
      if(debug){
        printAssembly("consensus2.fasta");
        BamPrinter bamPrinter = new BamPrinter();
        String bamPath = outPath + "mapped_reads2.bam";
        byte[] seq2 = getConsensusSimple(true, coverageThreshold, genome.seqB);
        Reference ref2 = new Reference(seq2, K, genome, false, threadNum);
        bamPrinter.printBAM(mappedData, ref2, bamPath);
      }
      logger.println("Simple consensus:");
      logger.println(new String(simpleConsensus));
      logger.println("2nd order consensus with filled gaps:");
      logger.println(new String(getConsensusSimple(true, coverageThreshold, genome.seqB)));
      logger.printf("Time iter 2, s: %d\n", (System.currentTimeMillis() - time)/1000);
      logger.println("Iteration 3");
      time = System.currentTimeMillis();
      //localize problem regions
      genome = new Reference(simpleConsensus, K, genome, false, threadNum);
      anchorSize = estimateAnchorSize(reads);
      ProblemIntervalBuilder problemIntervalBuilder = new ProblemIntervalBuilder(genome, mappedData,
              windowSize, identityThreshold, K, anchorSize, logger, threadNum);
      problemIntervals = problemIntervalBuilder.localizeProblemRegions();
      if(problemIntervals.length == 0){
        finalConsensus = new String(simpleConsensus);
        return finalConsensus;
      }
      // find all reads intersecting with edges of problem intervals
      selectGapIntersectingReads(reads);
      if(debug){
        printGaps();
      }
      //building problem reads index
      //mappedData = null;
      problemReadMapper.init(problemReads, mapper.counter, mapper.aligner, problemIntervals);
      problemReadMapper.mapProblemReadsLeft();
      while(true){
        //logger.println("Gap growing iteration " + iterNum);
        problemReadMapper.insertLoops();
        for(Interval interval : problemIntervals){
          if(interval.seqFinal == null && interval.finalizingReadsFound >= minFinalReads){
            mergeEdges(interval, mapper.aligner, identityThreshold,
                    minFinalReads, minIntersectLen);
          }
        }
        problemReadMapper.countLeftEdgeFreqs();
        // grow all problem region edges
        boolean leftUpdated = growLeftEdges();
        problemReadMapper.mapProblemReadsLeft();
        problemReadMapper.insertLoops();
        for(Interval interval : problemIntervals){
          if(interval.seqFinal == null && interval.finalizingReadsFound >= minFinalReads){
            mergeEdges(interval, mapper.aligner, identityThreshold,
                    minFinalReads, minIntersectLen);
          }
        }
        problemReadMapper.countRightEdgeFreqs();
        boolean rightUpdated = growRightEdges();
        problemReadMapper.mapProblemReadsRight();
        for(Interval interval : problemIntervals){
          if(interval.seqFinal == null && !interval.updatedLeft && !interval.updatedRight){
            mergeEdges(interval, mapper.aligner, identityThreshold,
                    minFinalReads, minIntersectLen);
          }
        }
        if(!leftUpdated && !rightUpdated){
          // no region edges were grown
          break;
        }
      }
      //finalizing all intervals
      for(Interval interval : problemIntervals){
        if(interval.seqFinal == null){
          mergeEdges(interval, mapper.aligner, identityThreshold,
                  minFinalReads, minIntersectLen);
        }
      }
      logger.println("Inserted:");
      for(Interval interval : problemIntervals){
        logger.println(interval.start + " " + interval.end + " " + interval.seqFinal);
      }
      // building final assembly
      genomeToAssembly = new int[genome.length];
      if(genomeIsFragmented){
        contigsFragmented = new ArrayList[genome.contigEnds.length];
        for(int i = 0; i < contigsFragmented.length; i++){
          contigsFragmented[i] = new ArrayList<>();
        }
        int fragmentID = 0;
        TByteArrayList builder = new TByteArrayList();
        int predPos = 0;
        int assemblyPos = 0;
        for(Interval interval : problemIntervals){
          while(interval.start >= genome.contigEnds[fragmentID]){
            int endPos = genome.contigEnds[fragmentID];
            if(builder.size() > 0){
              builder.add(simpleConsensus, predPos, endPos - predPos);
              contigsFragmented[fragmentID].add(builder.toArray());
              builder.clear();
            }else{
              contigsFragmented[fragmentID].add(Arrays.copyOfRange(simpleConsensus, predPos, endPos));
            }
            for(int i = predPos; i < endPos; i++){
              genomeToAssembly[i] = assemblyPos;
              assemblyPos ++;
            }
            predPos = endPos;
            fragmentID ++;
          }
          builder.add(simpleConsensus, predPos, interval.start - predPos);
          for(int i = predPos; i < interval.start; i++){
            genomeToAssembly[i] = assemblyPos;
            assemblyPos ++;
          }
          for(int i = interval.start; i < interval.end; i++){
            genomeToAssembly[i] = -1;
          }
          int i = interval.seqFinal.indexOf("------");
          if(i != -1){
            assemblyPos += interval.seqFinal.length() - 6;
            String left = interval.seqFinal.substring(0, i);
            String right = interval.seqFinal.substring(i + 6);
            builder.add(left.getBytes());
            if(builder.size() > 0){
              contigsFragmented[fragmentID].add(builder.toArray());
              builder.clear();
            }
            builder.add(right.getBytes());
            if(interval.end >= genome.contigEnds[fragmentID]){
              fragmentID ++;
            }
          }else{
            assemblyPos += interval.seqFinal.length();
            builder.add(interval.seqFinal.getBytes());
          }
          predPos = interval.end;
          while(fragmentID < genome.contigEnds.length && genome.contigEnds[fragmentID] <= predPos){
            fragmentID ++;
          }
        }
        if(fragmentID < genome.contigEnds.length){
          int endPos = genome.contigEnds[fragmentID];
          if(builder.size() > 0){
            builder.add(simpleConsensus, predPos, endPos - predPos);
            contigsFragmented[fragmentID].add(builder.toArray());
            builder.clear();
          }else{
            contigsFragmented[fragmentID].add(Arrays.copyOfRange(simpleConsensus, predPos, endPos));
          }
          for(int i = predPos; i < endPos; i++){
            genomeToAssembly[i] = assemblyPos;
            assemblyPos++;
          }
          predPos = endPos;
          fragmentID++;
          for(; fragmentID < genome.contigEnds.length; fragmentID++){
            endPos = genome.contigEnds[fragmentID];
            contigsFragmented[fragmentID].add(Arrays.copyOfRange(simpleConsensus, predPos, endPos));
            for(int i = predPos; i < endPos; i++){
              genomeToAssembly[i] = assemblyPos;
              assemblyPos++;
            }
            predPos = endPos;
          }
        }
        logger.println("Final consensus: ");
        StringBuilder fConsensus = new StringBuilder();
        TIntArrayList contigEnds = new TIntArrayList();
        fragmentEnds = new int[fragmentNames.length];
        for(int i = 0, end = 0; i < genome.fragmentNames.length; i++){
          logger.println(genome.fragmentNames[i]);
          for(int k = 0; k < contigsFragmented[i].size(); k++){
            byte[] subContig = contigsFragmented[i].get(k);
            String str = new String(subContig);
            fConsensus.append(str);
            end += subContig.length;
            contigEnds.add(end);
            logger.println(str);
          }
          fragmentEnds[i] = end;
        }
        this.contigEnds = contigEnds.toArray();
        finalConsensus = fConsensus.toString();
      }else{
        TByteArrayList builder = new TByteArrayList();
        int predPos = 0;
        int assemblyPos = 0;
        for(Interval interval : problemIntervals){
          builder.add(simpleConsensus, predPos, interval.start - predPos);
          for(int i = predPos; i < interval.start; i++){
            genomeToAssembly[i] = assemblyPos;
            assemblyPos ++;
          }
          for(int i = interval.start; i < interval.end; i++){
            genomeToAssembly[i] = -1;
          }
          if(interval.seqFinal.contains("------")){
            assemblyPos += interval.seqFinal.length() - 6;
          }else{
            assemblyPos += interval.seqFinal.length();
          }
          builder.add(interval.seqFinal.getBytes());
          predPos = interval.end;
        }
        builder.add(simpleConsensus, predPos, genome.length - predPos);
        for(int i = predPos; i < genome.length; i++){
          genomeToAssembly[i] = assemblyPos;
          assemblyPos ++;
        }
        finalConsensus = new String(builder.toArray());
        StringTokenizer tokenizer = new StringTokenizer(finalConsensus, "------");
        ArrayList<String> contigArr = new ArrayList<>();
        StringBuilder fConsensus = new StringBuilder();
        while(tokenizer.hasMoreTokens()){
          String token = tokenizer.nextToken();
          contigArr.add(token);
          fConsensus.append(token);
        }
        contigs = contigArr.toArray(new String[contigArr.size()]);
        finalConsensus = fConsensus.toString();
        logger.println("Final consensus: ");
        contigEnds = new int[contigs.length];
        int pos = 0;
        for(int i = 0; i < contigs.length; i++){
          String contig = contigs[i];
          logger.println(contig);
          pos += contig.length();
          contigEnds[i] = pos;
        }
      }
      logger.printf("Time iter 3, s: %d\n", (System.currentTimeMillis() - time)/1000);
    }catch(Exception e){
      e.printStackTrace();
    }
    return finalConsensus;
  }

  void printAssembly(String assemblyName) throws IOException{
    if(genomeIsFragmented){
      if(contigsFragmented == null){
        contigsFragmented = new ArrayList[fragmentNames.length];
        int prev = 0;
        for(int i = 0; i < fragmentEnds.length; i++){
          ArrayList<byte[]> contigs = new ArrayList<>();
          contigs.add(Arrays.copyOfRange(simpleConsensus, prev, fragmentEnds[i]));
          prev = fragmentEnds[i];
          contigsFragmented[i] = contigs;
        }
      }
      BufferedWriter writer = new BufferedWriter(new FileWriter(outPath + assemblyName));
      for(int i = 0; i < contigsFragmented.length; i++){
        if(contigsFragmented[i].size() == 1){
          byte[] subContig = contigsFragmented[i].get(0);
          String str = new String(subContig);
          writer.write(">" + fragmentNames[i].replace(' ', '_'));
          writer.newLine();
          writer.write(str);
          writer.newLine();
        }else{
          for(int k = 0; k < contigsFragmented[i].size(); k++){
            byte[] subContig = contigsFragmented[i].get(k);
            String str = new String(subContig);
            writer.write(">" + fragmentNames[i].replace(' ', '_') + "_" + Integer.toString(k + 1));
            writer.newLine();
            writer.write(str);
            writer.newLine();
          }
        }
      }
      writer.close();
    }else{
      if(contigs == null){
        contigs = new String[contigEnds.length];
        int prev = 0;
        for(int i = 0; i < contigEnds.length; i++){
          contigs[i] = new String(Arrays.copyOfRange(simpleConsensus, prev, contigEnds[i]));
          prev = contigEnds[i];
        }
      }
      BufferedWriter writer = new BufferedWriter(new FileWriter(outPath + assemblyName));
      for(int i = 0; i < contigs.length; i++){
        writer.write(">" + genomeName.replace(' ', '_') + "_" + Integer.toString(i + 1));
        writer.newLine();
        writer.write(contigs[i]);
        writer.newLine();
      }
      writer.close();
    }
  }

  private void adjustCoord(MappedRead mappedRead){
    mappedRead.end = genomeToAssembly[mappedRead.start + mappedRead.aln.end2 - 1] + 1;
    mappedRead.start = genomeToAssembly[mappedRead.start + mappedRead.aln.start2];
    mappedRead.aln.end2 -= mappedRead.aln.start2;
    mappedRead.aln.start2 = 0;
  }

  private void markRemapReads(){
    mappedData.needToRemap = new ArrayList<>();
    ArrayList<PairedRead> concordant = new ArrayList<>();
    ArrayList<PairedRead> leftMateMapped = new ArrayList<>();
    ArrayList<PairedRead> rightMateMapped = new ArrayList<>();
    ArrayList<MappedRead> mappedReads = new ArrayList<>();
    outer:
    for(PairedRead pairedRead: mappedData.concordant){
      for(Interval problemInterval: problemIntervals){
        if(pairedRead.r1.start + pairedRead.r1.aln.start2 < problemInterval.end &&
            pairedRead.r1.start + pairedRead.r1.aln.end2 > problemInterval.start){
          mappedData.needToRemap.add(pairedRead);
          continue outer;
        }
      }
      for(Interval problemInterval: problemIntervals){
        if(pairedRead.r2.start + pairedRead.r2.aln.start2 < problemInterval.end &&
            pairedRead.r2.start + pairedRead.r2.aln.end2 > problemInterval.start){
          mappedData.needToRemap.add(pairedRead);
          continue outer;
        }
      }
      concordant.add(pairedRead);
      adjustCoord(pairedRead.r1);
      adjustCoord(pairedRead.r2);
      mappedReads.add(pairedRead.r1);
      mappedReads.add(pairedRead.r2);
    }
    outer:
    for(PairedRead pairedRead: mappedData.leftMateMapped){
      if(!Arrays.equals(pairedRead.seq2, NULL_SEQ)){
        mappedData.needToRemap.add(pairedRead);
      }else{
        for(Interval problemInterval : problemIntervals){
          if(pairedRead.r1.start + pairedRead.r1.aln.start2 < problemInterval.end &&
              pairedRead.r1.start + pairedRead.r1.aln.end2 > problemInterval.start){
            mappedData.needToRemap.add(pairedRead);
            continue outer;
          }
        }
        leftMateMapped.add(pairedRead);
        adjustCoord(pairedRead.r1);
        mappedReads.add(pairedRead.r1);
      }
    }
    outer:
    for(PairedRead pairedRead: mappedData.rightMateMapped){
      if(!Arrays.equals(pairedRead.seq1, NULL_SEQ)){
        mappedData.needToRemap.add(pairedRead);
      }else{
        for(Interval problemInterval : problemIntervals){
          if(pairedRead.r2.start + pairedRead.r2.aln.start2 < problemInterval.end &&
              pairedRead.r2.start + pairedRead.r2.aln.end2 > problemInterval.start){
            mappedData.needToRemap.add(pairedRead);
            continue outer;
          }
        }
        rightMateMapped.add(pairedRead);
        adjustCoord(pairedRead.r2);
        mappedReads.add(pairedRead.r2);
      }
    }
    mappedData.needToRemap.addAll(mappedData.discordant);
    mappedData.needToRemap.addAll(mappedData.unmapped);
    mappedData.concordant = concordant;
    mappedData.leftMateMapped = leftMateMapped;
    mappedData.rightMateMapped = rightMateMapped;
    mappedData.mappedReads = mappedReads;
    mappedData.discordant.clear();
    mappedData.unmapped.clear();
  }

  void remapReadsToAssembly() throws InterruptedException{
    if(finalConsensus.equals("")){
      return;
    }
    Reference genome = new Reference(genomeName, finalConsensus, K, contigEnds, fragmentEnds,
        fragmentNames, genomeIsFragmented, threadNum);
    markRemapReads();
    mapper.mapReads(mappedData, genome);
  }

  public void assemble(Document jdomDocument) throws IOException, InterruptedException{
    long time = System.currentTimeMillis();
    DataReader dataReader = new DataReader(jdomDocument);
    Reference genome = new Reference(jdomDocument);
    ArrayList<PairedRead> reads = dataReader.readFilesWithReads();
    String finalConsensus = buildConsensus(genome, reads);
    if(finalConsensus.equals("")){
      printAssembly("assembly.fasta");
      logger.println("Total time: " + (System.currentTimeMillis() - time) / 1000);
      return;
    }
    genome = new Reference(genomeName, finalConsensus, K, contigEnds,
        fragmentEnds, fragmentNames, genomeIsFragmented, threadNum);
    remapReadsToAssembly();
    printAssembly("assembly.fasta");
    BamPrinter bamPrinter = new BamPrinter();
    String bamPath = outPath + "mapped_reads.bam";
    bamPrinter.printBAM(mappedData, genome, bamPath);
    //indexBam(bamPath);
    logger.println("Total time: " + (System.currentTimeMillis() - time) / 1000);
  }

  public static void main(String[] args){
    try{
      SAXBuilder jdomBuilder = new SAXBuilder();
      Document config = jdomBuilder.build(args[0]);
      ConsensusBuilderWithReassembling cBuilder = new ConsensusBuilderWithReassembling(config);
      cBuilder.assemble(config);
    }catch(Exception e){
      e.printStackTrace();
    }
  }

}
