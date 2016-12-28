import gnu.trove.list.array.TIntArrayList;
import org.jdom2.Document;
import org.jdom2.JDOMException;
import org.jdom2.input.SAXBuilder;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

/**
 * With postprocessing
 */
public class RefBasedAssembler extends Constants{

  private DataReader dataReader;
  private String outPath;
  private Document document;
  private Postprocessor postprocessor;
  private boolean useMajor;
  private ReferenceFinder refFinder;
  private Logger logger;
  private boolean usePostprocessing;

  public RefBasedAssembler(Document document){
    outPath = document.getRootElement().getChildText("OutPath");
    File out = new File(outPath);
    if(!out.exists()){
      out.mkdirs();
    }
    dataReader = new DataReader(document);
    usePostprocessing = Boolean.parseBoolean(
        document.getRootElement().getChild("Postprocessor").getChildText("Enabled"));
    postprocessor = new Postprocessor(document);
    useMajor = Boolean.parseBoolean(
        document.getRootElement().getChild("ReferenceSelector").getChildText("UseMajor"));
    refFinder = new ReferenceFinder(document);
    logger = Logger.getInstance(document);
    this.document = document;
  }

  private ConsensusBuilderWithReassembling buildUsingRef(Document document, Reference ref) throws IOException, InterruptedException{
    ConsensusBuilderWithReassembling consensusBuilder = new ConsensusBuilderWithReassembling(document);
    consensusBuilder.logger.println("Assembling using " + ref.name);
    ref.buildIndex(document);
    consensusBuilder.buildConsensus(ref, ref.reads);
    return consensusBuilder;
  }

  private void buildUsingRefNoPostproc(Document document, Reference ref) throws IOException, InterruptedException{
    ConsensusBuilderWithReassembling consensusBuilder = new ConsensusBuilderWithReassembling(document);
    ref.buildIndex(document);
    consensusBuilder.buildConsensus(ref, ref.reads);
    consensusBuilder.printAssembly(ref.name + "_assembly.fasta");
    Reference genome = new Reference(consensusBuilder);
    consensusBuilder.remapReadsToAssembly();
    BamPrinter bamPrinter = new BamPrinter();
    String bamPath = outPath + ref.name + "_mapped_reads.bam";
    bamPrinter.printBAM(consensusBuilder.mappedData, genome, bamPath);
  }

  private int[][] chooseIntervalsToMove(Interval[][] problemIntervals,
                                        ArrayList<Postprocessor.SendToInterval>[] intervals,
                                        int[][] coordConvert, String[] assemblies){
    int[][] problemIntervalsToMove = new int[intervals.length][];
    for(int i = 0; i < intervals.length; i++){
      int[] coordConv = coordConvert[i];
      problemIntervalsToMove[i] = new int[problemIntervals[i].length];
      outer:
      for(int j = 0; j < problemIntervals[i].length; j++){
        Interval interval = problemIntervals[i][j];
        int start;
        if(interval.start == 0){
          start = 0;
        }else{
          start = coordConv[interval.start - 1];
        }
        int end;
        if(interval.end == coordConv.length){
          end = assemblies[i].length();
        }else{
          end = coordConv[interval.end];
        }
        for(Postprocessor.SendToInterval sendToInterval: intervals[i]){
          if(start < sendToInterval.end && end > sendToInterval.start){
            problemIntervalsToMove[i][j] = sendToInterval.assemblyID;
//            System.out.printf("%d: %d - %d to %d\n", i, interval.start, interval.end,
//                sendToInterval.assemblyID);
            continue outer;
          }
        }
        problemIntervalsToMove[i][j] = -1;
      }
    }
    return problemIntervalsToMove;
  }

  private void adjustCoord(MappedRead mappedRead, int[] coordConvert){
    mappedRead.end = coordConvert[mappedRead.start + mappedRead.aln.end2 - 1] + 1;
    mappedRead.start = coordConvert[mappedRead.start + mappedRead.aln.start2];
    mappedRead.aln.end2 -= mappedRead.aln.start2;
    mappedRead.aln.start2 = 0;
  }

  private void resortConcordant(ArrayList<PairedRead> reads, Interval[] problemIntervals,
                                   int[] problemIntervalsToMove, int[] coordConv,
                                   MappedData mappedData_AfterCuts_i,
                                   ArrayList<Postprocessor.SendToInterval> sendToIntervals,
                                HashSet<String>[] mapped,
                                HashMap<String, PairedRead>[] needToRemap, int i){
    outer:
    for(PairedRead pairedRead: reads){
      // check if read is in problem interval
      boolean firstInProblemInterval = false;
      boolean secondInProblemInterval = false;
      for(int j = 0; j < problemIntervals.length; j++){
        Interval problemInterval = problemIntervals[j];
        if(pairedRead.r1.start + pairedRead.r1.aln.start2 < problemInterval.end &&
            pairedRead.r1.start + pairedRead.r1.aln.end2 > problemInterval.start){
          // it is
          int moveToID = problemIntervalsToMove[j];
          if(moveToID != -1){
            // if we need to move this interval to another assembly
            if(!mapped[moveToID].contains(pairedRead.name)){
              needToRemap[moveToID].put(pairedRead.name, pairedRead);
            }
            //MappedData mappedData1 = mappedDataAfterCuts[moveToID];
            //mappedData1.needToRemap.add(pairedRead);
            continue outer;
          }else{
            // we need to remap this read to current assembly
            firstInProblemInterval = true;
            //mappedData_AfterCuts_i.needToRemap.add(pairedRead);
          }
        }
        if(pairedRead.r2.start + pairedRead.r2.aln.start2 < problemInterval.end &&
            pairedRead.r2.start + pairedRead.r2.aln.end2 > problemInterval.start){
          // it is
          int moveToID = problemIntervalsToMove[j];
          if(moveToID != -1){
            // if we need to move this interval to another assembly
            if(!mapped[moveToID].contains(pairedRead.name)){
              needToRemap[moveToID].put(pairedRead.name, pairedRead);
            }
            //MappedData mappedData1 = mappedDataAfterCuts[moveToID];
            //mappedData1.needToRemap.add(pairedRead);
            continue outer;
          }else{
            // we need to remap this read to current assembly
            secondInProblemInterval = true;
            //mappedData_AfterCuts_i.needToRemap.add(pairedRead);
          }

        }
      }
      // this read pair is not in problem interval -> compute it's new coord on curr assembly
      int start = coordConv[pairedRead.r1.start + pairedRead.r1.aln.start2];
      int end = coordConv[pairedRead.r1.start + pairedRead.r1.aln.end2 - 1] + 1;
      // check if this read should be moved to new assembly
      for(Postprocessor.SendToInterval sendToInterval: sendToIntervals){
        if(start < sendToInterval.end && end > sendToInterval.start){
//            MappedData mappedDataTo = mappedDataAfterCuts[sendToInterval.assemblyID];
//            mappedDataTo.needToRemap.add(pairedRead);
          if(!mapped[sendToInterval.assemblyID].contains(pairedRead.name)){
            needToRemap[sendToInterval.assemblyID].put(pairedRead.name, pairedRead);
          }
          continue outer;
        }
      }
      start = coordConv[pairedRead.r2.start + pairedRead.r2.aln.start2];
      end = coordConv[pairedRead.r2.start + pairedRead.r2.aln.end2 - 1] + 1;
      // check if this read should be moved to new assembly
      for(Postprocessor.SendToInterval sendToInterval: sendToIntervals){
        if(start < sendToInterval.end && end > sendToInterval.start){
//            MappedData mappedDataTo = mappedDataAfterCuts[sendToInterval.assemblyID];
//            mappedDataTo.needToRemap.add(pairedRead);
          if(!mapped[sendToInterval.assemblyID].contains(pairedRead.name)){
            needToRemap[sendToInterval.assemblyID].put(pairedRead.name, pairedRead);
          }
          continue outer;
        }
      }
      if(firstInProblemInterval || secondInProblemInterval){
        needToRemap[i].put(pairedRead.name, pairedRead);
      }else{
        // this read stay here and not remapped
        adjustCoord(pairedRead.r1, coordConv);
        adjustCoord(pairedRead.r2, coordConv);
        mappedData_AfterCuts_i.concordant.add(pairedRead);
        mappedData_AfterCuts_i.mappedReads.add(pairedRead.r1);
        mappedData_AfterCuts_i.mappedReads.add(pairedRead.r2);
        mapped[i].add(pairedRead.name);
        needToRemap[i].remove(pairedRead.name);
      }
    }
  }

  private void resortDiscordant(ArrayList<PairedRead> reads, Interval[] problemIntervals,
                                int[] problemIntervalsToMove, int[] coordConv,
                                ArrayList<Postprocessor.SendToInterval> sendToIntervals,
                                HashSet<String>[] mapped,
                                HashMap<String, PairedRead>[] needToRemap, int i){
    outer:
    for(PairedRead pairedRead: reads){
      // check if read is in problem interval
      for(int j = 0; j < problemIntervals.length; j++){
        Interval problemInterval = problemIntervals[j];
        if(pairedRead.r1.start + pairedRead.r1.aln.start2 < problemInterval.end &&
            pairedRead.r1.start + pairedRead.r1.aln.end2 > problemInterval.start){
          // it is
          int moveToID = problemIntervalsToMove[j];
          if(moveToID != -1){
            // if we need to move this interval to another assembly
//              MappedData mappedData1 = mappedDataAfterCuts[moveToID];
//              mappedData1.needToRemap.add(pairedRead);
            if(!mapped[moveToID].contains(pairedRead.name)){
              needToRemap[moveToID].put(pairedRead.name, pairedRead);
            }
            continue outer;
          }else{
            // we need to remap this read to current assembly
            //mappedData_AfterCuts_i.needToRemap.add(pairedRead);
          }
        }
        if(pairedRead.r2.start + pairedRead.r2.aln.start2 < problemInterval.end &&
            pairedRead.r2.start + pairedRead.r2.aln.end2 > problemInterval.start){
          // it is
          int moveToID = problemIntervalsToMove[j];
          if(moveToID != -1){
            // if we need to move this interval to another assembly
//              MappedData mappedData1 = mappedDataAfterCuts[moveToID];
//              mappedData1.needToRemap.add(pairedRead);
            if(!mapped[moveToID].contains(pairedRead.name)){
              needToRemap[moveToID].put(pairedRead.name, pairedRead);
            }
            continue outer;
          }else{
            // we need to remap this read to current assembly
            //mappedData_AfterCuts_i.needToRemap.add(pairedRead);
          }
        }
      }
      // this read pair is no in problem interval -> compute it's new coord on curr assembly
      int start = coordConv[pairedRead.r1.start + pairedRead.r1.aln.start2];
      int end = coordConv[pairedRead.r1.start + pairedRead.r1.aln.end2 - 1] + 1;
      // check if this read should be moved to new assembly
      for(Postprocessor.SendToInterval sendToInterval: sendToIntervals){
        if(start < sendToInterval.end && end > sendToInterval.start){
//            MappedData mappedDataTo = mappedDataAfterCuts[sendToInterval.assemblyID];
//            mappedDataTo.needToRemap.add(pairedRead);
          if(!mapped[sendToInterval.assemblyID].contains(pairedRead.name)){
            needToRemap[sendToInterval.assemblyID].put(pairedRead.name, pairedRead);
          }
          continue outer;
        }
      }
      start = coordConv[pairedRead.r2.start + pairedRead.r2.aln.start2];
      end = coordConv[pairedRead.r2.start + pairedRead.r2.aln.end2 - 1] + 1;
      // check if this read should be moved to new assembly
      for(Postprocessor.SendToInterval sendToInterval: sendToIntervals){
        if(start < sendToInterval.end && end > sendToInterval.start){
//            MappedData mappedDataTo = mappedDataAfterCuts[sendToInterval.assemblyID];
//            mappedDataTo.needToRemap.add(pairedRead);
          if(!mapped[sendToInterval.assemblyID].contains(pairedRead.name)){
            needToRemap[sendToInterval.assemblyID].put(pairedRead.name, pairedRead);
          }
          continue outer;
        }
      }
      // we remap it anyway
      //mappedData_AfterCuts_i.needToRemap.add(pairedRead);
      needToRemap[i].put(pairedRead.name, pairedRead);
    }
  }

  private void resortLeft(ArrayList<PairedRead> reads, Interval[] problemIntervals,
                                int[] problemIntervalsToMove, int[] coordConv,
                                MappedData mappedData_AfterCuts_i,
                                ArrayList<Postprocessor.SendToInterval> sendToIntervals,
                                HashSet<String>[] mapped,
                                HashMap<String, PairedRead>[] needToRemap, int i){
    outer:
    for(PairedRead pairedRead: reads){
      // check if read is in problem interval
      boolean firstInProblemInterval = false;
      for(int j = 0; j < problemIntervals.length; j++){
        Interval problemInterval = problemIntervals[j];
        if(pairedRead.r1.start + pairedRead.r1.aln.start2 < problemInterval.end &&
            pairedRead.r1.start + pairedRead.r1.aln.end2 > problemInterval.start){
          // it is
          int moveToID = problemIntervalsToMove[j];
          if(moveToID != -1){
            // if we need to move this interval to another assembly
//              MappedData mappedData1 = mappedDataAfterCuts[moveToID];
//              mappedData1.needToRemap.add(pairedRead);
            if(!mapped[moveToID].contains(pairedRead.name)){
              needToRemap[moveToID].put(pairedRead.name, pairedRead);
            }
            continue outer;
          }else{
            // we need to remap this read to current assembly
            firstInProblemInterval = true;
            //mappedData_AfterCuts_i.needToRemap.add(pairedRead);
          }
        }
      }
      // this read pair is no in problem interval -> compute it's new coord on curr assembly
      int start = coordConv[pairedRead.r1.start + pairedRead.r1.aln.start2];
      int end = coordConv[pairedRead.r1.start + pairedRead.r1.aln.end2 - 1] + 1;
      // check if this read should be moved to new assembly
      for(Postprocessor.SendToInterval sendToInterval: sendToIntervals){
        if(start < sendToInterval.end && end > sendToInterval.start){
//            MappedData mappedDataTo = mappedDataAfterCuts[sendToInterval.assemblyID];
//            mappedDataTo.needToRemap.add(pairedRead);
          if(!mapped[sendToInterval.assemblyID].contains(pairedRead.name)){
            needToRemap[sendToInterval.assemblyID].put(pairedRead.name, pairedRead);
          }
          continue outer;
        }
      }
      if(!Arrays.equals(pairedRead.seq2, NULL_SEQ)){
        // pair is present -> remap
        //mappedData_AfterCuts_i.needToRemap.add(pairedRead);
        needToRemap[i].put(pairedRead.name, pairedRead);
      }else{
        // one read from pair present
        if(firstInProblemInterval){
          //mappedData_AfterCuts_i.needToRemap.add(pairedRead);
          needToRemap[i].put(pairedRead.name, pairedRead);
        }else{
          mappedData_AfterCuts_i.leftMateMapped.add(pairedRead);
          mappedData_AfterCuts_i.mappedReads.add(pairedRead.r1);
          adjustCoord(pairedRead.r1, coordConv);
          mapped[i].add(pairedRead.name);
          needToRemap[i].remove(pairedRead.name);
        }
      }
    }
  }

  private void resortRight(ArrayList<PairedRead> reads, Interval[] problemIntervals,
                          int[] problemIntervalsToMove, int[] coordConv,
                          MappedData mappedData_AfterCuts_i,
                          ArrayList<Postprocessor.SendToInterval> sendToIntervals,
                          HashSet<String>[] mapped,
                          HashMap<String, PairedRead>[] needToRemap, int i){
    outer:
    for(PairedRead pairedRead: reads){
      boolean secondInProblemInterval = false;
      // check if read is in problem interval
      for(int j = 0; j < problemIntervals.length; j++){
        Interval problemInterval = problemIntervals[j];
        if(pairedRead.r2.start + pairedRead.r2.aln.start2 < problemInterval.end &&
            pairedRead.r2.start + pairedRead.r2.aln.end2 > problemInterval.start){
          // it is
          int moveToID = problemIntervalsToMove[j];
          if(moveToID != -1){
            // if we need to move this interval to another assembly
//              MappedData mappedData1 = mappedDataAfterCuts[moveToID];
//              mappedData1.needToRemap.add(pairedRead);
            if(!mapped[moveToID].contains(pairedRead.name)){
              needToRemap[moveToID].put(pairedRead.name, pairedRead);
            }
            continue outer;
          }else{
            // we need to remap this read to current assembly
            secondInProblemInterval = true;
            //mappedData_AfterCuts_i.needToRemap.add(pairedRead);
          }
        }
      }
      // this read pair is no in problem interval -> compute it's new coord on curr assembly
      int start = coordConv[pairedRead.r2.start + pairedRead.r2.aln.start2];
      int end = coordConv[pairedRead.r2.start + pairedRead.r2.aln.end2 - 1] + 1;
      // check if this read should be moved to new assembly
      for(Postprocessor.SendToInterval sendToInterval: sendToIntervals){
        if(start < sendToInterval.end && end > sendToInterval.start){
//            MappedData mappedDataTo = mappedDataAfterCuts[sendToInterval.assemblyID];
//            mappedDataTo.needToRemap.add(pairedRead);
          if(!mapped[sendToInterval.assemblyID].contains(pairedRead.name)){
            needToRemap[sendToInterval.assemblyID].put(pairedRead.name, pairedRead);
          }
          continue outer;
        }
      }
      if(!Arrays.equals(pairedRead.seq1, NULL_SEQ)){
        // pair is present -> remap
        //mappedData_AfterCuts_i.needToRemap.add(pairedRead);
        needToRemap[i].put(pairedRead.name, pairedRead);
      }else{
        // one read from pair present
        if(secondInProblemInterval){
          //mappedData_AfterCuts_i.needToRemap.add(pairedRead);
          needToRemap[i].put(pairedRead.name, pairedRead);
        }else{
          mappedData_AfterCuts_i.rightMateMapped.add(pairedRead);
          mappedData_AfterCuts_i.mappedReads.add(pairedRead.r2);
          adjustCoord(pairedRead.r2, coordConv);
          mapped[i].add(pairedRead.name);
          needToRemap[i].remove(pairedRead.name);
        }
      }
    }
  }

  private MappedData[] resortReads(Interval[][] problemIntervals,
                                  ArrayList<Postprocessor.SendToInterval>[] intervals,
                                  MappedData[] mappedData, int[][] coordConvert,
                                  int[][] problemIntervalsToMove){
    MappedData[] mappedDataAfterCuts = new MappedData[intervals.length];
    HashMap<String, PairedRead>[] needToRemap = new HashMap[intervals.length];
    HashSet<String>[] mapped = new HashSet[intervals.length];
    for(int i = 0; i < mappedDataAfterCuts.length; i++){
      mappedDataAfterCuts[i] = new MappedData();
      needToRemap[i] = new HashMap<>();
      mapped[i] = new HashSet<>();
    }
    for(int i = 0; i < intervals.length; i++){
      ArrayList<Postprocessor.SendToInterval> sendToIntervals = intervals[i];
      MappedData mappedData_i = mappedData[i];
      MappedData mappedData_AfterCuts_i = mappedDataAfterCuts[i];
      int[] coordConv = coordConvert[i];
      resortConcordant(mappedData_i.concordant, problemIntervals[i], problemIntervalsToMove[i], coordConv,
          mappedData_AfterCuts_i, sendToIntervals, mapped, needToRemap, i);
      resortDiscordant(mappedData_i.discordant, problemIntervals[i], problemIntervalsToMove[i], coordConv,
          sendToIntervals, mapped, needToRemap, i);
      resortLeft(mappedData_i.leftMateMapped, problemIntervals[i], problemIntervalsToMove[i], coordConv,
          mappedData_AfterCuts_i, sendToIntervals, mapped, needToRemap, i);
      resortRight(mappedData_i.rightMateMapped, problemIntervals[i], problemIntervalsToMove[i], coordConv,
          mappedData_AfterCuts_i, sendToIntervals, mapped, needToRemap, i);
      for(PairedRead pairedRead: mappedData_i.unmapped){
        needToRemap[i].put(pairedRead.name, pairedRead);
      }
    }
    for(int i = 0; i < intervals.length; i++){
      mappedDataAfterCuts[i].needToRemap.addAll(needToRemap[i].values());
    }
    return mappedDataAfterCuts;
  }

  private void reformatAssembly(ConsensusBuilderWithReassembling consensusBuilder,
                                ArrayList<Postprocessor.SendToInterval> sendToIntervals){
    byte[] finalConsensus = consensusBuilder.finalConsensus.getBytes();
    int[] coordTransform = new int[finalConsensus.length];
    if(consensusBuilder.genomeIsFragmented){
      ArrayList<byte[]>[] contigsFragmented = new ArrayList[consensusBuilder.contigEnds.length];
      for(int i = 0; i < contigsFragmented.length; i++){
        contigsFragmented[i] = new ArrayList<>();
      }
      int[] fragmentEnds = new int[consensusBuilder.contigsFragmented.length];
      for(int i = 0; i < consensusBuilder.contigsFragmented.length; i++){
        for(byte[] contig: consensusBuilder.contigsFragmented[i]){
          fragmentEnds[i] += contig.length;
        }
      }
      int fragmentID = 0;
      int contigID = 0;
      int predPos = 0;
      int assemblyPos = 0;
      for(Postprocessor.SendToInterval interval : sendToIntervals){
        while(interval.start >= consensusBuilder.contigEnds[contigID]){
          int endPos = consensusBuilder.contigEnds[contigID];
          while(endPos > fragmentEnds[fragmentID]){
            fragmentID ++;
          }
          if(endPos > predPos){
            contigsFragmented[fragmentID].add(Arrays.copyOfRange(finalConsensus, predPos, endPos));
          }
          for(int i = predPos; i < endPos; i++){
            coordTransform[i] = assemblyPos;
            assemblyPos ++;
          }
          predPos = endPos;
          contigID ++;
        }
        contigsFragmented[fragmentID].add(Arrays.copyOfRange(finalConsensus, predPos, interval.start));
        for(int i = predPos; i < interval.start; i++){
          coordTransform[i] = assemblyPos;
          assemblyPos ++;
        }
        for(int i = interval.start; i < interval.end; i++){
          coordTransform[i] = -1;
        }
        predPos = interval.end;
        while(consensusBuilder.contigEnds[contigID] < predPos){
          contigID ++;
        }
      }
      int endPos = consensusBuilder.contigEnds[contigID];
      while(endPos > fragmentEnds[fragmentID]){
        fragmentID ++;
      }
      if(endPos > predPos){
        contigsFragmented[fragmentID].add(Arrays.copyOfRange(finalConsensus, predPos, endPos));
      }
      for(int i = predPos; i < endPos; i++){
        coordTransform[i] = assemblyPos;
        assemblyPos ++;
      }
      predPos = endPos;
      contigID ++;
      for(;contigID < consensusBuilder.contigEnds.length; contigID ++){
        endPos = consensusBuilder.contigEnds[contigID];
        while(endPos > fragmentEnds[fragmentID]){
          fragmentID ++;
        }
        contigsFragmented[fragmentID].add(Arrays.copyOfRange(finalConsensus, predPos, endPos));
        for(int i = predPos; i < endPos; i++){
          coordTransform[i] = assemblyPos;
          assemblyPos ++;
        }
        predPos = endPos;
      }
      consensusBuilder.logger.println("Final consensus after postprocessing: ");
      StringBuilder fConsensus = new StringBuilder();
      TIntArrayList contigEnds = new TIntArrayList();
      for(int i = 0, j = 0; i < consensusBuilder.contigEnds.length; i++){
        consensusBuilder.logger.println(consensusBuilder.fragmentNames[i]);
        for(int k = 0; k < contigsFragmented[i].size(); k++){
          byte[] subContig = contigsFragmented[i].get(k);
          String str = new String(subContig);
          fConsensus.append(str);
          j += subContig.length;
          contigEnds.add(j);
          consensusBuilder.logger.println(str);
        }
      }
      consensusBuilder.contigsFragmented = contigsFragmented;
      consensusBuilder.finalConsensus = fConsensus.toString();
      consensusBuilder.contigEnds = contigEnds.toArray();
      consensusBuilder.genomeToAssembly = coordTransform;
    }else{
      int predPos = 0;
      ArrayList<byte[]> contigs = new ArrayList<>();
      int contigID = 0;
      int assemblyPos = 0;
      for(Postprocessor.SendToInterval interval : sendToIntervals){
        while(interval.start >= consensusBuilder.contigEnds[contigID]){
          int endPos = consensusBuilder.contigEnds[contigID];
          if(endPos > predPos){
            contigs.add(Arrays.copyOfRange(finalConsensus, predPos, endPos));
          }
          for(int i = predPos; i < endPos; i++){
            coordTransform[i] = assemblyPos;
            assemblyPos ++;
          }
          predPos = endPos;
          contigID ++;
        }
        contigs.add(Arrays.copyOfRange(finalConsensus, predPos, interval.start));
        for(int i = predPos; i < interval.start; i++){
          coordTransform[i] = assemblyPos;
          assemblyPos ++;
        }
        for(int i = interval.start; i < interval.end; i++){
          coordTransform[i] = -1;
        }
        predPos = interval.end;
        while(consensusBuilder.contigEnds[contigID] < predPos){
          contigID ++;
        }
      }
      int endPos = consensusBuilder.contigEnds[contigID];
      if(endPos > predPos){
        contigs.add(Arrays.copyOfRange(finalConsensus, predPos, endPos));
      }
      for(int i = predPos; i < endPos; i++){
        coordTransform[i] = assemblyPos;
        assemblyPos ++;
      }
      predPos = endPos;
      contigID ++;
      for(;contigID < consensusBuilder.contigEnds.length; contigID ++){
        endPos = consensusBuilder.contigEnds[contigID];
        contigs.add(Arrays.copyOfRange(finalConsensus, predPos, endPos));
        for(int i = predPos; i < endPos; i++){
          coordTransform[i] = assemblyPos;
          assemblyPos ++;
        }
        predPos = endPos;
      }
      consensusBuilder.contigs = new String[contigs.size()];
      consensusBuilder.logger.println("Final consensus after postprocessing: ");
      StringBuilder fConsensus = new StringBuilder();
      TIntArrayList contigEnds = new TIntArrayList();
      for(int i = 0, j = 0; i < contigs.size(); i++){
        byte[] contig = contigs.get(i);
        String s = new String(contig);
        consensusBuilder.contigs[i] = s;
        fConsensus.append(s);
        j += contig.length;
        contigEnds.add(j);
        consensusBuilder.logger.println(s);
      }
      consensusBuilder.finalConsensus = fConsensus.toString();
      consensusBuilder.contigEnds = contigEnds.toArray();
      consensusBuilder.genomeToAssembly = coordTransform;
    }
  }

  public void assemble() throws IOException, InterruptedException, InvocationTargetException, InstantiationException, IllegalAccessException, JDOMException{
    long time = System.currentTimeMillis();
    ArrayList<PairedRead> reads = dataReader.readFilesWithReads();
    ArrayList<Reference> selectedRefs = refFinder.selectReferences(reads);
    if(selectedRefs.size() == 0){
      Reference genome = new Reference(document);
      genome.reads = reads;
      selectedRefs.add(genome);
    }
    if(useMajor){
      int maxReadNum = 0;
      Reference maxRef = null;
      for(Reference seq : selectedRefs){
        if(seq.reads.size() > maxReadNum){
          maxReadNum = seq.reads.size();
          maxRef = seq;
        }
      }
      BufferedWriter writer = new BufferedWriter(new FileWriter(outPath +
          "selected_reference.fasta"));
      writer.write(">" + maxRef.name + "\n");
      writer.write(maxRef.seq + "\n");
      writer.close();
      selectedRefs.clear();
      selectedRefs.add(maxRef);
    }
    if(usePostprocessing){
      String[] assemblies = new String[selectedRefs.size()];
      Reference[] selRefs = selectedRefs.toArray(new Reference[assemblies.length]);
      MappedData[] mappedData = new MappedData[assemblies.length];
      int[][] coordConvert = new int[assemblies.length][];
      Interval[][] problemIntervals = new Interval[assemblies.length][];
      ConsensusBuilderWithReassembling[] cBuilders = new ConsensusBuilderWithReassembling[assemblies.length];
      int[] readNums = new int[assemblies.length];
      for(int i = 0; i < selRefs.length; i++){
        ConsensusBuilderWithReassembling consensusBuilderGC = buildUsingRef(document, selRefs[i]);
        assemblies[i] = consensusBuilderGC.finalConsensus;
        mappedData[i] = consensusBuilderGC.mappedData;
        coordConvert[i] = consensusBuilderGC.genomeToAssembly;
        problemIntervals[i] = consensusBuilderGC.problemIntervals;
        cBuilders[i] = consensusBuilderGC;
        readNums[i] = consensusBuilderGC.reads.size();
        Reference ref = selRefs[i];
        consensusBuilderGC.printAssembly(ref.name.replace(' ', '_') + "_assembly_no_postproc.fasta");
      }
      ArrayList<Postprocessor.SendToInterval>[] intervals =
          postprocessor.findRegionsToCut(assemblies, selRefs, readNums);
      int[][] problemIntervalsToMove = chooseIntervalsToMove(problemIntervals, intervals,
          coordConvert, assemblies);
      MappedData[] mappedDataAfterCuts = resortReads(problemIntervals, intervals, mappedData,
          coordConvert, problemIntervalsToMove);
      for(int i = 0; i < assemblies.length; i++){
        ConsensusBuilderWithReassembling cBuilder = cBuilders[i];
        reformatAssembly(cBuilder, intervals[i]);
        cBuilder.mappedData = mappedDataAfterCuts[i];
        for(MappedRead mappedRead : cBuilder.mappedData.mappedReads){
          adjustCoord(mappedRead, cBuilder.genomeToAssembly);
        }
        Reference ref = selRefs[i];
        Reference genome = new Reference(cBuilder);
        cBuilder.mapper.mapReads(cBuilder.mappedData, genome);
        recomputeAssemblyConsensus(cBuilder);
        cBuilder.printAssembly(ref.name.replace(' ', '_') + "_assembly.fasta");
        BamPrinter bamPrinter = new BamPrinter();
        String bamPath = outPath + ref.name.replace(' ', '_') + "_mapped_reads.bam";
        bamPrinter.printBAM(mappedDataAfterCuts[i], genome, bamPath);
        //cBuilder.indexBam(bamPath);
      }
    }else{
      for(Reference ref : selectedRefs){
        buildUsingRefNoPostproc(document, ref);
      }
    }
    logger.println("Total time: " + (System.currentTimeMillis() - time)/1000);
  }

  private void recomputeAssemblyConsensus(ConsensusBuilderWithReassembling consensusBuilder) throws InterruptedException{
//    Reference genome = new Reference(consensusBuilder.finalConsensus);
    byte[] consensus = consensusBuilder.getConsensusSimple(true, 1, consensusBuilder.finalConsensus.getBytes());
    consensusBuilder.finalConsensus = new String(consensus);
    if(consensusBuilder.genomeIsFragmented){
      ArrayList<byte[]>[] contigsFragmented =
          new ArrayList[consensusBuilder.contigsFragmented.length];
      for(int i = 0, k = 0, predPos = 0; i < contigsFragmented.length; i++){
        ArrayList<byte[]> list = new ArrayList<>();
        ArrayList<byte[]> old = consensusBuilder.contigsFragmented[i];
        for(int j = 0; j < old.size(); j++, k++){
          int contigEnd = consensusBuilder.contigEnds[k];
          list.add(consensusBuilder.finalConsensus.substring(predPos, contigEnd).getBytes());
          predPos = contigEnd;
        }
        contigsFragmented[i] = list;
      }
      consensusBuilder.contigsFragmented = contigsFragmented;
    }else{
      String[] contigs = new String[consensusBuilder.contigs.length];
      for(int i = 0, predPos = 0; i < consensusBuilder.contigEnds.length; i++){
        int contigEnd = consensusBuilder.contigEnds[i];
        contigs[i] = consensusBuilder.finalConsensus.substring(predPos, contigEnd);
        predPos = contigEnd;
      }
      consensusBuilder.contigs = contigs;
    }
  }

  private static void printUsage(){
    System.out.println("Usage: java -cp ./ViRGA.jar RefBasedAssembler pathToConfigFile");
  }

  public static void main(String[] args){
    try{
      SAXBuilder jdomBuilder = new SAXBuilder();
      if(args.length == 0 || args[0].equals("-h")){
        printUsage();
        return;
      }
      if(args.length != 1){
        System.out.println("Wrong parameter number!");
        printUsage();
        return;
      }
      Document config = jdomBuilder.build(args[0]);
      boolean selectRefs = Boolean.parseBoolean(
          config.getRootElement().getChild("ReferenceSelector").getChildText("Enabled"));
      if(selectRefs){
        RefBasedAssembler assembler = new RefBasedAssembler(config);
        assembler.assemble();
      }else{
        ConsensusBuilderWithReassembling cBuilder = new ConsensusBuilderWithReassembling(config);
        cBuilder.assemble(config);
      }
    }catch(Exception e){
      e.printStackTrace();
    }
  }

}
