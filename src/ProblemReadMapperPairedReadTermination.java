import gnu.trove.iterator.TObjectIntIterator;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TObjectIntHashMap;
import org.jdom2.Document;
import org.jdom2.Element;

import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Created by Геннадий on 08.12.2014.
 */
class ProblemReadMapperPairedReadTermination extends ProblemReadMapper {

  private LinkedList<ProblemRead>[] workingLists;
  private float identityThreshold;
  private int coverageThreshold;
  private float gapMergeThreshold;
  private int minLen;
  private int K;

  ProblemReadMapperPairedReadTermination(Document document) {
    super(document);
    Element element = document.getRootElement().getChild("ConsensusBuilder").getChild("Reassembler");
    identityThreshold = Float.parseFloat(element.getChildText("IdentityThreshold"));
    coverageThreshold = Integer.parseInt(element.getChildText("MinTerminatingSequenceCoverage"));
    gapMergeThreshold = Float.parseFloat(element.getChildText("PairReadTerminationThreshold"));
    minLen = Integer.parseInt(element.getChildText("MinReadLength"));
    K = KMerCounter.getInstance(document).K;
  }

  private MappedRead mapReadVeryFast(KMerCounter counter, byte[] read, int genLen,
                                     HashMap<String, int[]> index) {
    MappedRead res = counter.mapReadToRegion(read, index, 0, genLen);
    if (res == null) {
      res = new MappedRead();
      res.count = 0;
    }
    res.seq = read;
    return res;
  }

  private void balanceWorkingLists() {
    int total = 0;
    for (LinkedList<ProblemRead> list : workingLists) {
      total += list.size();
    }
    int av = total / workingLists.length;
    for (int i = 0; i < workingLists.length; i++) {
      LinkedList<ProblemRead> list_i = workingLists[i];
      int excess = list_i.size() - av;
      if (excess > 0) {
        for (int j = 0; j < i; j++) {
          LinkedList<ProblemRead> list_j = workingLists[j];
          int lack = av - list_j.size();
          if (lack > 0) {
            Iterator<ProblemRead> iter = list_i.iterator();
            if (excess > lack) {
              for (int k = 0; k < lack; k++) {
                ProblemRead read = iter.next();
                list_j.add(read);
                iter.remove();
              }
              excess -= lack;
            } else {
              for (int k = 0; k < excess; k++) {
                ProblemRead read = iter.next();
                list_j.add(read);
                iter.remove();
              }
              break;
            }
          }
        }
      }
      if (excess > 0) {
        for (int j = i + 1; j < workingLists.length; j++) {
          LinkedList<ProblemRead> list_j = workingLists[j];
          int lack = av - list_j.size();
          if (lack > 0) {
            Iterator<ProblemRead> iter = list_i.iterator();
            if (excess > lack) {
              for (int k = 0; k < lack; k++) {
                ProblemRead read = iter.next();
                list_j.add(read);
                iter.remove();
              }
              excess -= lack;
            } else {
              for (int k = 0; k < excess; k++) {
                ProblemRead read = iter.next();
                list_j.add(read);
                iter.remove();
              }
              break;
            }
          }
        }
      }
    }
  }

  private class ReadMappingLeft implements Runnable {

    private int k;
    private KMerCounter counter;
    private Aligner aligner;
    int mappedReads = 0;
    int mappedForward = 0;
    int mappedReverse = 0;
    HashSet<String>[] leftReads;
    HashSet<String>[] rightReads;

    ReadMappingLeft(KMerCounter counter, Aligner aligner,
                           int k) {
      this.counter = counter;
      this.aligner = aligner;
      this.k = k;
      if(debug){
        leftReads = new HashSet[problemIntervals.length];
        rightReads = new HashSet[problemIntervals.length];
        for(int i = 0; i < problemIntervals.length; i++){
          leftReads[i] = new HashSet<String>();
          rightReads[i] = new HashSet<String>();
        }
      }
    }

    @Override
    public void run() {
      Iterator<ProblemRead> iter = workingLists[k].iterator();
      while (iter.hasNext()) {
        ProblemRead read = iter.next();
        short gapID = read.gapID;
        ProblemInterval gap = problemIntervals[gapID];
        if(gap.merge || gap.seqFinal != null) {//gap.seqFinal != null
          iter.remove();
          continue;
        }
        short sizeConcat = (short) gap.concat.length;
        if (read.count == -1) {
          read.setMapping(mapReadVeryFast(counter, read.seq, sizeConcat, gap.index));
          read.lastLen = sizeConcat;
        }
        boolean mapsToConcat = read.count > read.cutoff;
        if (!mapsToConcat &&
            sizeConcat - read.lastLen > read.cutoff - read.count) {
          read.setMapping(mapReadVeryFast(counter, read.seq, sizeConcat, gap.index));
          read.lastLen = sizeConcat;
          mapsToConcat = read.count > read.cutoff;
        }
        if (mapsToConcat) {
          mappedReads++;
          if (read.reverse == 0) {
            mappedForward++;
          } else {
            mappedReverse++;
          }
          if (read.aln == null) {
            read.countLeft = (short) counter.computeKMerCount(read.seq, gap.indexLeft, 0, gap.left.length);
            read.lastLenLeft = (short) gap.left.length;
            read.countRight = (short) counter.computeKMerCount(read.seq, gap.indexRight, 0, gap.right.length);
            read.lastLenRight = (short) gap.right.length;
            read.alignSW(gap, aligner);
            if (read.aln.junctionLeft != -1 && read.aln.junctionRight != -1) {
              if ((float) read.aln.leftIdentity / read.aln.leftAlnLen <= identityThreshold ||
                  read.aln.leftIdentity < read.seq.length / 3) {
                read.aln.junctionLeft = -1;
                read.aln.leftAlnLen = 0;
              }
              if ((float) read.aln.rightIdentity / read.aln.rightAlnLen <= identityThreshold ||
                  read.aln.rightIdentity < read.seq.length / 3) {
                read.aln.junctionRight = -1;
                read.aln.rightAlnLen = 0;
              }
            }
            if (read.aln.junctionLeft == -1) {
              if (read.aln.junctionRight == -1) {
                read.mappingType = -1;
                read.aln = null;
              } else {
                read.mappingType = 1;
                short sizeRight = (short) gap.right.length;
                MappedRead mRead = mapReadVeryFast(counter, read.seq, sizeRight, gap.indexRight);
                read.countRight = (short) mRead.count;
                read.lastLenRight = sizeRight;
                read.setMapping(mRead);
                read.alignSWRight(gap, aligner);
              }
            } else {
              if (read.aln.junctionRight == -1) {
                read.mappingType = 0;
                short sizeLeft = (short) gap.left.length;
                MappedRead mRead = mapReadVeryFast(counter, read.seq, sizeLeft, gap.indexLeft);
                read.countLeft = (short) mRead.count;
                read.lastLenLeft = sizeLeft;
                read.setMapping(mRead);
                read.alignSWLeft(gap, aligner);
              } else {
                read.mappingType = 2;
              }
            }
          } else {
            if (gap.updatedLeft) {
              if (read.mappingType == 0) {
                if (read.aln.junctionLeft != -1) {
                  read.extendAlignmentLeft(gap);
                } else {
                  short sizeLeft = (short) gap.left.length;
                  MappedRead mRead = mapReadVeryFast(counter, read.seq, sizeLeft, gap.indexLeft);
                  if (mRead.count > read.countLeft) {
                    read.countLeft = (short) mRead.count;
                    read.lastLenLeft = sizeLeft;
                    read.setMapping(mRead);
                    read.alignSWLeft(gap, aligner);
                  }
                }
              } else if (read.mappingType == 2) {
                read.extendAlignmentLeft(gap);
              } else if (read.mappingType == 1) {
                if (read.aln.junctionRight != -1 && read.aln.junctionRight > read.seq.length / 3) {
                  short sizeLeft = (short) gap.left.length;
                  MappedRead mRead = mapReadVeryFast(counter, read.seq, sizeLeft, gap.indexLeft);
                  if (mRead.count > read.countLeft) {
                    read.countLeft = (short) mRead.count;
                    read.lastLenLeft = sizeLeft;
                    ProblemRead pRead = new ProblemRead(read.name, read.n, read.seq, gapID, read.reverse);
                    pRead.setMapping(mapReadVeryFast(counter, read.seq, sizeConcat, gap.index));
                    pRead.alignSW(gap, aligner);
                    if (pRead.aln.junctionLeft != -1 && pRead.aln.junctionRight != -1) {
                      if ((float) pRead.aln.leftIdentity / pRead.aln.leftAlnLen <= identityThreshold ||
                          pRead.aln.leftIdentity < pRead.seq.length / 3) {
                        pRead.aln.junctionLeft = -1;
                        pRead.aln.leftAlnLen = 0;
                      }
                      if ((float) pRead.aln.rightIdentity / pRead.aln.rightAlnLen <= identityThreshold ||
                          pRead.aln.rightIdentity < pRead.seq.length / 3) {
                        pRead.aln.junctionRight = -1;
                        pRead.aln.rightAlnLen = 0;
                      }
                    }
                    if (pRead.aln.junctionLeft != -1 && pRead.aln.junctionRight != -1) {
                      read.start = pRead.start;
                      read.aln = pRead.aln;
                      read.end = pRead.end;
                      read.mappingType = 2;
                    }
                  }
                }
              }
            }
          }
          if (read.mappingType == 0) {
            if (read.seq.length == read.aln.end1 ||
                read.end + read.seq.length - read.aln.end1 <= gap.left.length) {
              iter.remove();
            }
          } else if (read.mappingType == 1) {
            if (read.aln.start1 == 0 || read.start - read.aln.start1 > 0) {
              iter.remove();
            }
          }
          if(debug){
            if(read.mappingType == 0){
              rightReads[gapID].add(read.name + "_" + read.n);
            }else if(read.mappingType == 2){
              leftReads[gapID].add(read.name + "_" + read.n);
              rightReads[gapID].add(read.name + "_" + read.n);
            }
          }
        }
      }
    }

  }

  private class ReadMappingRight implements Runnable {

    private int k;
    private KMerCounter counter;
    private Aligner aligner;
    int mappedReads = 0;
    int mappedForward = 0;
    int mappedReverse = 0;
    HashSet<String>[] leftReads;
    HashSet<String>[] rightReads;

    ReadMappingRight(KMerCounter counter, Aligner aligner,
                            int k) {
      this.counter = counter;
      this.aligner = aligner;
      this.k = k;
      if(debug){
        leftReads = new HashSet[problemIntervals.length];
        rightReads = new HashSet[problemIntervals.length];
        for(int i = 0; i < problemIntervals.length; i++){
          leftReads[i] = new HashSet<String>();
          rightReads[i] = new HashSet<String>();
        }
      }
    }

    @Override
    public void run() {
      Iterator<ProblemRead> iter = workingLists[k].iterator();
      while (iter.hasNext()) {
        ProblemRead read = iter.next();
        short gapID = read.gapID;
        ProblemInterval gap = problemIntervals[gapID];
        if (gap.merge || gap.seqFinal != null) {//gap.seqFinal != null
          iter.remove();
          continue;
        }
        short sizeConcat = (short) gap.concat.length;
        if (read.count == -1) {
          read.setMapping(mapReadVeryFast(counter, read.seq, sizeConcat, gap.index));
          read.lastLen = sizeConcat;
        }
        boolean mapsToConcat = read.count > read.cutoff;
        if (!mapsToConcat &&
            sizeConcat - read.lastLen > read.cutoff - read.count) {
          read.setMapping(mapReadVeryFast(counter, read.seq, sizeConcat, gap.index));
          read.lastLen = sizeConcat;
          mapsToConcat = read.count > read.cutoff;
        }
        if (mapsToConcat) {
          mappedReads++;
          if (read.reverse == 0) {
            mappedForward++;
          } else {
            mappedReverse++;
          }
          if (read.aln == null) {
            read.countLeft = (short) counter.computeKMerCount(read.seq, gap.indexLeft, 0, gap.left.length);
            read.lastLenLeft = (short) gap.left.length;
            read.countRight = (short) counter.computeKMerCount(read.seq, gap.indexRight, 0, gap.right.length);
            read.lastLenRight = (short) gap.right.length;
            read.alignSW(gap, aligner);
            if (read.aln.junctionLeft != -1 && read.aln.junctionRight != -1) {
              if ((float) read.aln.leftIdentity / read.aln.leftAlnLen <= identityThreshold ||
                  read.aln.leftIdentity < read.seq.length / 3) {
                read.aln.junctionLeft = -1;
                read.aln.leftAlnLen = 0;
              }
              if ((float) read.aln.rightIdentity / read.aln.rightAlnLen <= identityThreshold ||
                  read.aln.rightIdentity < read.seq.length / 3) {
                read.aln.junctionRight = -1;
                read.aln.rightAlnLen = 0;
              }
            }
            if (read.aln.junctionLeft == -1) {
              if (read.aln.junctionRight == -1) {
                read.mappingType = -1;
              } else {
                read.mappingType = 1;
                short sizeRight = (short) gap.right.length;
                MappedRead mRead = mapReadVeryFast(counter, read.seq, sizeRight, gap.indexRight);
                read.countRight = (short) mRead.count;
                read.lastLenRight = sizeRight;
                read.setMapping(mRead);
                read.alignSWRight(gap, aligner);
              }
            } else {
              if (read.aln.junctionRight == -1) {
                read.mappingType = 0;
                short sizeLeft = (short) gap.left.length;
                MappedRead mRead = mapReadVeryFast(counter, read.seq, sizeLeft, gap.indexLeft);
                read.countLeft = (short) mRead.count;
                read.lastLenLeft = sizeLeft;
                read.setMapping(mRead);
                read.alignSWLeft(gap, aligner);
              } else {
                read.mappingType = 2;
              }
            }
          }else{
            if (gap.updatedRight) {
              if (read.mappingType == 1) {
                if (read.aln.junctionRight != -1) {
                  read.extendAlignmentRight(gap);
                } else {
                  MappedRead mRead = mapReadVeryFast(counter, read.seq, sizeConcat, gap.indexRight);
                  if (mRead.count > read.countRight) {
                    read.countRight = (short) mRead.count;
                    read.lastLenRight = (short) gap.right.length;
                    read.setMapping(mRead);
                    read.alignSWRight(gap, aligner);
                  }
                }
              } else if (read.mappingType == 2) {
                read.extendAlignmentRight(gap);
              } else if (read.mappingType == 0) {
                if (read.aln.junctionLeft != -1 &&
                    read.seq.length - read.aln.junctionLeft > read.seq.length / 3) {
                  MappedRead mRead = mapReadVeryFast(counter, read.seq, sizeConcat, gap.indexRight);
                  if (mRead.count > read.countRight) {
                    read.countRight = (short) mRead.count;
                    read.lastLenRight = (short) gap.right.length;
                    ProblemRead pRead = new ProblemRead(read.name, read.n, read.seq, gapID, read.reverse);
                    pRead.setMapping(mapReadVeryFast(counter, read.seq, sizeConcat, gap.index));
                    pRead.alignSW(gap, aligner);
                    if (pRead.aln.junctionLeft != -1 && pRead.aln.junctionRight != -1) {
                      if ((float) pRead.aln.leftIdentity / pRead.aln.leftAlnLen <= identityThreshold ||
                          pRead.aln.leftIdentity < pRead.seq.length / 3) {
                        pRead.aln.junctionLeft = -1;
                        pRead.aln.leftAlnLen = 0;
                      }
                      if ((float) pRead.aln.rightIdentity / pRead.aln.rightAlnLen <= identityThreshold ||
                          pRead.aln.rightIdentity < pRead.seq.length / 3) {
                        pRead.aln.junctionRight = -1;
                        pRead.aln.rightAlnLen = 0;
                      }
                    }
                    if (pRead.aln.junctionRight != -1 && pRead.aln.junctionLeft != -1) {
                      read.start = pRead.start;
                      read.aln = pRead.aln;
                      read.end = pRead.end;
                      read.mappingType = 2;
                    }
                  }
                }
              }
            }
          }
          if (read.mappingType == 0) {
            if (read.seq.length == read.aln.end1 ||
                read.end + read.seq.length - read.aln.end1 <= gap.left.length) {
              iter.remove();
            }
          } else if (read.mappingType == 1) {
            if (read.aln.start1 == 0 || read.start - read.aln.start1 > 0) {
              iter.remove();
            }
          }
          if(debug){
            if(read.mappingType == 1){
              rightReads[gapID].add(read.name + "_" + read.n);
            }else if(read.mappingType == 2){
              leftReads[gapID].add(read.name + "_" + read.n);
              rightReads[gapID].add(read.name + "_" + read.n);
            }
          }
        }
      }
    }

  }

  @Override
  void init(ArrayList<ProblemRead> problemReads, KMerCounter counter, Aligner aligner,
                      ProblemInterval[] problemIntervals) {
    super.init(problemReads, counter, aligner, problemIntervals);
    workingLists = new LinkedList[threadNum];
    for (int i = 0; i < threadNum; i++) {
      workingLists[i] = new LinkedList<>();
    }
    int k = 0;
    for (ProblemRead read : problemReads) {
      if (read.gapID != -1) {
        workingLists[k % threadNum].add(read);
        k++;
      }
    }
  }

  void mapProblemReadsLeft()
      throws InterruptedException{
    int totalReads = 0;
    for (LinkedList<ProblemRead> list : workingLists) {
      totalReads += list.size();
    }
    int mappedReads = 0;
    int mappedForward = 0;
    int mappedReverse = 0;
    ReadMappingLeft[] tasks = new ReadMappingLeft[threadNum];
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    for (int i = 0; i < threadNum; i++){
      ReadMappingLeft task = new ReadMappingLeft(counter, aligner, i);
      executor.execute(task);
      tasks[i] = task;
    }
    executor.shutdown();
    executor.awaitTermination(100, TimeUnit.HOURS);
    for (int n = 0; n < threadNum; n++){
      ReadMappingLeft task = tasks[n];
      mappedReads += task.mappedReads;
      mappedForward += task.mappedForward;
      mappedReverse += task.mappedReverse;
      if(debug){
        for(int i = 0; i < problemIntervals.length; i++){
          leftReads[i].addAll(task.leftReads[i]);
          rightReads[i].addAll(task.rightReads[i]);
        }
      }
    }
    balanceWorkingLists();
    if(debug){
      logger.printf("Total problem reads: %d\n", totalReads);
      logger.printf("Total problem reads mapped left: %1.2f, forward: %1.2f, reverse: %1.2f\n",
          (float) mappedReads/totalReads, (float) mappedForward/mappedReads,
          (float) mappedReverse/mappedReads);
    }
  }

  void mapProblemReadsRight()
      throws InterruptedException{
    int totalReads = 0;
    for (LinkedList<ProblemRead> list : workingLists){
      totalReads += list.size();
    }
    int mappedReads = 0;
    int mappedForward = 0;
    int mappedReverse = 0;
    ReadMappingRight[] tasks = new ReadMappingRight[threadNum];
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    for (int i = 0; i < threadNum; i++) {
      ReadMappingRight task = new ReadMappingRight(counter, aligner, i);
      executor.execute(task);
      tasks[i] = task;
    }
    executor.shutdown();
    executor.awaitTermination(100, TimeUnit.HOURS);
    for(int n = 0; n < threadNum; n++){
      ReadMappingRight task = tasks[n];
      mappedReads += task.mappedReads;
      mappedForward += task.mappedForward;
      mappedReverse += task.mappedReverse;
      if(debug){
        for(int i = 0; i < problemIntervals.length; i++){
          leftReads[i].addAll(task.leftReads[i]);
          rightReads[i].addAll(task.rightReads[i]);
        }
      }
    }
    balanceWorkingLists();
    if(debug){
      logger.printf("Total problem reads: %d\n", totalReads);
      logger.printf("Total problem reads mapped right: %1.2f, forward: %1.2f, reverse: %1.2f\n",
          (float) mappedReads/totalReads, (float) mappedForward/mappedReads,
          (float) mappedReverse/mappedReads);
    }
  }

  private class LeftEdgeFreqCounting implements Runnable {

    private int k;
    public int[][] gapLeftCounts;

    public LeftEdgeFreqCounting(int k) {
      this.k = k;
    }

    @Override
    public void run() {
      gapLeftCounts = new int[problemIntervals.length][4];
      for (ProblemRead read : workingLists[k]) {
        short gapID = read.gapID;
        if (!problemIntervals[gapID].updateCount) {
          continue;
        }
        boolean mapsToConcat = read.count >= read.cutoff;
        if (mapsToConcat) {
          if (read.mappingType == 0) {
            if (read.aln.junctionLeft == -1) {
              continue;
            }
            if ((float) read.aln.leftIdentity / read.aln.leftAlnLen > identityThreshold) {
              byte c = read.seq[read.aln.junctionLeft + 1];
              if (c != 'N') {
                gapLeftCounts[gapID][nToI[c]]++;
              }
            }
          }
        }
      }
    }

  }

  @Override
  void countLeftEdgeFreqs() {
    try {
      LeftEdgeFreqCounting[] tasks = new LeftEdgeFreqCounting[threadNum];
      ExecutorService executor = Executors.newFixedThreadPool(threadNum);
      for (int i = 0; i < threadNum; i++) {
        LeftEdgeFreqCounting task = new LeftEdgeFreqCounting(i);
        executor.execute(task);
        tasks[i] = task;
      }
      executor.shutdown();
      executor.awaitTermination(100, TimeUnit.HOURS);
      for (int n = 0; n < threadNum; n++) {
        LeftEdgeFreqCounting task = tasks[n];
        for (int i = 0; i < problemIntervals.length; i++) {
          ProblemInterval gap = problemIntervals[i];
          if (!gap.updateCount) {
            continue;
          }
          for (int j = 0; j < 4; j++) {
            gap.leftCounts[j] += task.gapLeftCounts[i][j];
          }
        }
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  private class RightEdgeFreqCounting implements Runnable {

    private int k;
    public int[][] gapRightCounts;

    public RightEdgeFreqCounting(int k) {
      this.k = k;
    }

    @Override
    public void run() {
      gapRightCounts = new int[problemIntervals.length][4];
      for (ProblemRead read : workingLists[k]) {
        short gapID = read.gapID;
        if (!problemIntervals[gapID].updateCount) {
          continue;
        }
        boolean mapsToConcat = read.count >= read.cutoff;
        if (mapsToConcat) {
          if (read.mappingType == 1) {
            if (read.aln.junctionRight == -1) {
              continue;
            }
            if ((float) read.aln.rightIdentity / read.aln.rightAlnLen > identityThreshold) {
              byte c = read.seq[read.aln.junctionRight - 1];
              if (c != 'N') {
                gapRightCounts[gapID][nToI[c]]++;
              }
            }
          }
        }
      }
    }

  }

  @Override
  void countRightEdgeFreqs() {
    try {
      RightEdgeFreqCounting[] tasks = new RightEdgeFreqCounting[threadNum];
      ExecutorService executor = Executors.newFixedThreadPool(threadNum);
      for (int i = 0; i < threadNum; i++) {
        RightEdgeFreqCounting task = new RightEdgeFreqCounting(i);
        executor.execute(task);
        tasks[i] = task;
      }
      executor.shutdown();
      executor.awaitTermination(100, TimeUnit.HOURS);
      for (int n = 0; n < threadNum; n++) {
        RightEdgeFreqCounting task = tasks[n];
        for (int i = 0; i < problemIntervals.length; i++) {
          ProblemInterval gap = problemIntervals[i];
          if (!gap.updateCount) {
            continue;
          }
          for (int j = 0; j < 4; j++) {
            gap.rightCounts[j] += task.gapRightCounts[i][j];
          }
        }
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }


  private class InsertsCountingTwoSide implements Runnable {

    private int k;
    public TObjectIntHashMap<String>[] insertCounts;

    public InsertsCountingTwoSide(int k) {
      this.k = k;
    }

    @Override
    public void run() {
      insertCounts = new TObjectIntHashMap[problemIntervals.length];
      for (int i = 0; i < problemIntervals.length; i++) {
        insertCounts[i] = new TObjectIntHashMap<>();
      }
      for (ProblemRead read : workingLists[k]) {
        short gapID = read.gapID;
        boolean mapsToConcat = read.count >= read.cutoff;
        if (mapsToConcat) {
          if (read.mappingType == 2) {
            boolean leftIdentity = (float) read.aln.leftIdentity / read.aln.leftAlnLen > identityThreshold;
            boolean rightIdentity = (float) read.aln.rightIdentity / read.aln.rightAlnLen > identityThreshold;
            if (read.aln.leftIdentity >= read.seq.length / 3 &&
                read.aln.rightIdentity >= read.seq.length / 3 && leftIdentity && rightIdentity) {
              TObjectIntHashMap<String> insertToCount = insertCounts[gapID];
              String insert = new String(Arrays.copyOfRange(read.seq, read.aln.junctionLeft + 1, read.aln.junctionRight));
              int count = insertToCount.get(insert);
              insertToCount.put(insert, count + 1);
            }
          }
        }
      }
    }

  }

  private TObjectIntHashMap<String>[] countInsertsTwoSide()
      throws InterruptedException {
    InsertsCountingTwoSide[] tasks = new InsertsCountingTwoSide[threadNum];
    //ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    Thread[] threads = new Thread[threadNum];
    for (int i = 0; i < threadNum; i++) {
      InsertsCountingTwoSide task = new InsertsCountingTwoSide(i);
      //executor.execute(task);
      threads[i] = new Thread(task);
      threads[i].start();
      tasks[i] = task;
    }
    boolean allComplete = false;
    while (!allComplete) {
      allComplete = true;
      for (int i = 0; i < threadNum; i++) {
        if (threads[i].isAlive()) {
          allComplete = false;
        }
      }
      Thread.currentThread().sleep(1);
    }
    //executor.shutdown();
    //executor.awaitTermination(100, TimeUnit.HOURS);
    TObjectIntHashMap<String>[] res = tasks[0].insertCounts;
    for (int n = 1; n < threadNum; n++) {
      InsertsCountingTwoSide task = tasks[n];
      for (int i = 0; i < problemIntervals.length; i++) {
        TObjectIntHashMap<String> insertToCounts = task.insertCounts[i];
        TObjectIntIterator<String> iter = insertToCounts.iterator();
        TObjectIntHashMap<String> res_i = res[i];
        while (iter.hasNext()) {
          iter.advance();
          res_i.adjustOrPutValue(iter.key(), iter.value(), iter.value());
        }
      }
    }
    return res;
  }

  private class PairCounter implements Runnable{

    int gapID;
    int totalActive;
    int oneUnmapped;
    int k;
    String insert;
    ArrayList<ProblemRead> unmapped;

    PairCounter(int gapID, int k, String insert){
      this.gapID = gapID;
      this.k = k;
      this.insert = insert;
    }

    HashMap<String, int[]> indexConcat(String seq){
      HashMap<String, TIntArrayList> index = new HashMap<>();
      for(int i = 0; i <= seq.length() - K; i++){
        String s = seq.substring(i, i + K);
        TIntArrayList l = index.get(s);
        if(l == null){
          l = new TIntArrayList();
          index.put(s, l);
        }
        l.add(i);
      }
      HashMap<String, int[]> map = new HashMap<>();
      for(Map.Entry<String, TIntArrayList> entry: index.entrySet()){
        map.put(entry.getKey(), entry.getValue().toArray());
      }
      return map;
    }

    @Override
    public void run(){
      unmapped = new ArrayList<>();
      ProblemInterval gap = problemIntervals[gapID];
      HashMap<String, int[]> index = indexConcat(new String(gap.left) + insert + new String(gap.right));
      int sizeConcat = gap.concat.length + insert.length();
      for(ProblemRead read: workingLists[k]){
        if(read.gapID != gapID || read.mate == null ||
            read.mate.seq.length < minLen){
          continue;
        }
        if(read.mappingType != -1){
          totalActive ++;
          continue;
        }
        if(sizeConcat - read.lastLen >= read.cutoff - read.count){
          int newCount = counter.computeKMerCount(read.seq, index, 0, sizeConcat);
          if(newCount > read.cutoff){
            totalActive ++;
            continue;
          }
        }
        if(read.mate instanceof ProblemRead){
          ProblemRead mate = (ProblemRead)read.mate;
          if(mate.gapID == gapID){
            if(mate.mappingType == -1){
              if(sizeConcat - mate.lastLen >= mate.cutoff - mate.count){
                int newCount = counter.computeKMerCount(mate.seq, index, 0, sizeConcat);
                if(newCount > mate.cutoff){
                  oneUnmapped ++;
                  totalActive ++;
                }
              }
            }else{
              oneUnmapped ++;
              totalActive ++;
            }
          }
        }else{
          oneUnmapped ++;
          totalActive ++;
        }
      }
    }

  }

  private boolean isMostOfPairsMapped(int gapID, String insert, int threadNum) throws InterruptedException{
    PairCounter[] tasks = new PairCounter[threadNum];
    ExecutorService executor = Executors.newFixedThreadPool(threadNum);
    for(int i = 0; i < threadNum; i++){
      PairCounter task = new PairCounter(gapID, i, insert);
      executor.execute(task);
      tasks[i] = task;
    }
    executor.shutdown();
    executor.awaitTermination(100, TimeUnit.HOURS);
    int totalActive = 0;
    int oneUnmapped = 0;
    ProblemInterval gap = problemIntervals[gapID];
    try{
      for(PairCounter task : tasks){
        totalActive += task.totalActive;
        oneUnmapped += task.oneUnmapped;
      }
    }catch(Exception e){
      e.printStackTrace();
    }
    if(debug){
      logger.println("Gap " + gap.start + " " + gap.end + ": total = " + totalActive +
          " one = " + oneUnmapped);
    }
    return (float)oneUnmapped/totalActive < gapMergeThreshold;
  }

  void insertLoops(){
    try{
      TObjectIntHashMap<String>[] inserts = countInsertsTwoSide();
      for(int i = 0; i < problemIntervals.length; i++){
        ProblemInterval gap = problemIntervals[i];
        gap.updateCount = false;
        if(gap.seqFinal == null){
          TObjectIntHashMap inserts_i = inserts[i];
          String insert = null;
          int maxCount = 0;
          TObjectIntIterator<String> iter = inserts_i.iterator();
          while(iter.hasNext()){
            iter.advance();
            if(iter.value() > maxCount){
              maxCount = iter.value();
              insert = iter.key();
            }
          }
          if(maxCount > coverageThreshold && isMostOfPairsMapped(i, insert, threadNum)){
            byte[] right = new byte[gap.right.length + insert.length()];
            System.arraycopy(gap.right, 0, right, insert.length(), gap.right.length);
            System.arraycopy(insert.getBytes(), 0, right, 0, insert.length());
            gap.right = right;
            byte[] concat = new byte[gap.concat.length + insert.length()];
            System.arraycopy(gap.left, 0, concat, 0, gap.left.length);
            System.arraycopy(gap.right, 0, concat, gap.junction + 1, gap.right.length);
            gap.concat = concat;
            gap.finalizingReadsFound = maxCount;
          }else{
            gap.updateCount = true;
          }
        }
      }
    }catch(Exception e){
      e.printStackTrace();
    }
  }




}
