import gnu.trove.list.array.TByteArrayList;
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
 * Created by Gennady on 19.05.2017.
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
    if(refFinder.debug && ref.problemReads.size() != 0){
      MappedData garbage = consensusBuilder.mapper.mapReads(ref.problemReads, genome);
      bamPath = outPath + ref.name + "_garbage.bam";
      bamPrinter.printBAM(garbage, genome, bamPath);
    }
  }


  private void adjustCoord(MappedRead mappedRead, int[] coordConvert){
    mappedRead.end = coordConvert[mappedRead.start + mappedRead.aln.end2 - 1] + 1;
    mappedRead.start = coordConvert[mappedRead.start + mappedRead.aln.start2];
    mappedRead.aln.end2 -= mappedRead.aln.start2;
    mappedRead.aln.start2 = 0;
  }

  private int computeIntersectionLen(int start, int end, RepeatInterval repeatInterval){
    int intersectionLen = 0;
    if(start < repeatInterval.start){
      if(end < repeatInterval.end){
        intersectionLen = end - repeatInterval.start;
      }else{
        intersectionLen = repeatInterval.end - repeatInterval.start;
      }
    }else{
      if(end < repeatInterval.end){
        intersectionLen = end - start;
      }else{
        intersectionLen = repeatInterval.end - start;
      }
    }
    return intersectionLen;
  }

  private void resortConcordant(ArrayList<PairedRead> reads,
                                MappedData mappedData_AfterCuts_i,
                                ArrayList<RepeatInterval> repeatIntervals,
                                HashSet<String>[] mapped,
                                HashMap<String, PairedRead>[] needToRemap, int i){
    outer:
    for(PairedRead pairedRead: reads){
      // compute new coord on curr assembly
      int start = pairedRead.r1.start + pairedRead.r1.aln.start2;
      int end = pairedRead.r1.start + pairedRead.r1.aln.end2;
      // check if this read should be moved to new assembly
      for(RepeatInterval repeatInterval : repeatIntervals){
        if(start < repeatInterval.end && end > repeatInterval.start){
          if(!mapped[repeatInterval.assemblyID].contains(pairedRead.name) &&
              computeIntersectionLen(start, end, repeatInterval) > 0.5*pairedRead.r1.aln.length){
            needToRemap[repeatInterval.assemblyID].put(pairedRead.name, pairedRead);
          }else{
            needToRemap[i].put(pairedRead.name, pairedRead);
          }
          continue outer;
        }
      }
      start = pairedRead.r2.start + pairedRead.r2.aln.start2;
      end = pairedRead.r2.start + pairedRead.r2.aln.end2;
      // check if this read should be moved to new assembly
      for(RepeatInterval repeatInterval : repeatIntervals){
        if(start < repeatInterval.end && end > repeatInterval.start){
          if(!mapped[repeatInterval.assemblyID].contains(pairedRead.name) &&
              computeIntersectionLen(start, end, repeatInterval) > 0.5*pairedRead.r2.aln.length){
            needToRemap[repeatInterval.assemblyID].put(pairedRead.name, pairedRead);
          }else{
            needToRemap[i].put(pairedRead.name, pairedRead);
          }
          continue outer;
        }
      }
      // this read stay here and not remapped
      mappedData_AfterCuts_i.concordant.add(pairedRead);
      mappedData_AfterCuts_i.mappedReads.add(pairedRead.r1);
      mappedData_AfterCuts_i.mappedReads.add(pairedRead.r2);
    }
  }

  private void resortDiscordant(ArrayList<PairedRead> reads,
                                ArrayList<RepeatInterval> repeatIntervals,
                                HashSet<String>[] mapped,
                                HashMap<String, PairedRead>[] needToRemap, int i){
    outer:
    for(PairedRead pairedRead: reads){
      // this read pair is no in problem interval -> compute it's new coord on curr assembly
      int start = pairedRead.r1.start + pairedRead.r1.aln.start2;
      int end = pairedRead.r1.start + pairedRead.r1.aln.end2;
      // check if this read should be moved to new assembly
      for(RepeatInterval repeatInterval : repeatIntervals){
        if(start < repeatInterval.end && end > repeatInterval.start){
          if(!mapped[repeatInterval.assemblyID].contains(pairedRead.name) &&
              computeIntersectionLen(start, end, repeatInterval) > 0.5*pairedRead.r1.aln.length){
            needToRemap[repeatInterval.assemblyID].put(pairedRead.name, pairedRead);
          }else{
            needToRemap[i].put(pairedRead.name, pairedRead);
          }
          continue outer;
        }
      }
      start = pairedRead.r2.start + pairedRead.r2.aln.start2;
      end = pairedRead.r2.start + pairedRead.r2.aln.end2;
      // check if this read should be moved to new assembly
      for(RepeatInterval repeatInterval : repeatIntervals){
        if(start < repeatInterval.end && end > repeatInterval.start){
          if(!mapped[repeatInterval.assemblyID].contains(pairedRead.name) &&
              computeIntersectionLen(start, end, repeatInterval) > 0.5*pairedRead.r2.aln.length){
            needToRemap[repeatInterval.assemblyID].put(pairedRead.name, pairedRead);
          }else{
            needToRemap[i].put(pairedRead.name, pairedRead);
          }
          continue outer;
        }
      }
      // we remap it anyway
      needToRemap[i].put(pairedRead.name, pairedRead);
    }
  }

  private void resortLeft(ArrayList<PairedRead> reads,
                          MappedData mappedData_AfterCuts_i,
                          ArrayList<RepeatInterval> repeatIntervals,
                          HashSet<String>[] mapped,
                          HashMap<String, PairedRead>[] needToRemap, int i){
    outer:
    for(PairedRead pairedRead: reads){
      // this read pair is no in problem interval -> compute it's new coord on curr assembly
      int start = pairedRead.r1.start + pairedRead.r1.aln.start2;
      int end = pairedRead.r1.start + pairedRead.r1.aln.end2;
      // check if this read should be moved to new assembly
      for(RepeatInterval repeatInterval : repeatIntervals){
        if(start < repeatInterval.end && end > repeatInterval.start){
          if(!mapped[repeatInterval.assemblyID].contains(pairedRead.name) &&
              computeIntersectionLen(start, end, repeatInterval) > 0.5*pairedRead.r1.aln.length){
            needToRemap[repeatInterval.assemblyID].put(pairedRead.name, pairedRead);
          }else{
            needToRemap[i].put(pairedRead.name, pairedRead);
          }
          continue outer;
        }
      }
      if(!Arrays.equals(pairedRead.seq2, NULL_SEQ)){
        // pair is present -> remap
        needToRemap[i].put(pairedRead.name, pairedRead);
      }else{
        // one read from pair present
        mappedData_AfterCuts_i.leftMateMapped.add(pairedRead);
        mappedData_AfterCuts_i.mappedReads.add(pairedRead.r1);
      }
    }
  }

  private void resortRight(ArrayList<PairedRead> reads,
                           MappedData mappedData_AfterCuts_i,
                           ArrayList<RepeatInterval> repeatIntervals,
                           HashSet<String>[] mapped,
                           HashMap<String, PairedRead>[] needToRemap, int i){
    outer:
    for(PairedRead pairedRead: reads){
      // this read pair is no in problem interval -> compute it's new coord on curr assembly
      int start = pairedRead.r2.start + pairedRead.r2.aln.start2;
      int end = pairedRead.r2.start + pairedRead.r2.aln.end2;
      // check if this read should be moved to new assembly
      for(RepeatInterval repeatInterval : repeatIntervals){
        if(start < repeatInterval.end && end > repeatInterval.start){
//            MappedData mappedDataTo = mappedDataAfterCuts[repeatInterval.assemblyID];
//            mappedDataTo.needToRemap.add(pairedRead);
          if(!mapped[repeatInterval.assemblyID].contains(pairedRead.name) &&
              computeIntersectionLen(start, end, repeatInterval) > 0.5*pairedRead.r2.aln.length){
            needToRemap[repeatInterval.assemblyID].put(pairedRead.name, pairedRead);
          }else{
            needToRemap[i].put(pairedRead.name, pairedRead);
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
        mappedData_AfterCuts_i.rightMateMapped.add(pairedRead);
        mappedData_AfterCuts_i.mappedReads.add(pairedRead.r2);
      }
    }
  }

  private MappedData[] resortReads(ArrayList<RepeatInterval>[] intervals,
                                   MappedData[] mappedData){
    MappedData[] mappedDataAfterCuts = new MappedData[intervals.length];
    HashMap<String, PairedRead>[] needToRemap = new HashMap[intervals.length];
    HashSet<String>[] mapped = new HashSet[intervals.length];
    for(int i = 0; i < mappedDataAfterCuts.length; i++){
      mappedDataAfterCuts[i] = new MappedData();
      needToRemap[i] = new HashMap<>();
      HashSet<String> mapped_i = new HashSet<>();
      for(MappedRead mappedRead: mappedData[i].mappedReads){
        mapped_i.add(mappedRead.name);
      }
      mapped[i] = mapped_i;
    }
    for(int i = 0; i < intervals.length; i++){
      ArrayList<RepeatInterval> repeatIntervals = intervals[i];
      MappedData mappedData_i = mappedData[i];
      MappedData mappedData_AfterCuts_i = mappedDataAfterCuts[i];
      resortConcordant(mappedData_i.concordant,
          mappedData_AfterCuts_i, repeatIntervals, mapped, needToRemap, i);
      resortDiscordant(mappedData_i.discordant,
          repeatIntervals, mapped, needToRemap, i);
      resortLeft(mappedData_i.leftMateMapped,
          mappedData_AfterCuts_i, repeatIntervals, mapped, needToRemap, i);
      resortRight(mappedData_i.rightMateMapped,
          mappedData_AfterCuts_i, repeatIntervals, mapped, needToRemap, i);
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
                                ArrayList<RepeatInterval> repeatIntervals){
    byte[] finalConsensus = consensusBuilder.finalConsensus.getBytes();
    if(finalConsensus.length == 0){
      return;
    }
    int[] coordTransform = new int[finalConsensus.length];
    if(consensusBuilder.genomeIsFragmented){
      ArrayList<byte[]>[] contigsFragmented = new ArrayList[consensusBuilder.contigEnds.length];
      for(int i = 0; i < contigsFragmented.length; i++){
        contigsFragmented[i] = new ArrayList<>();
      }
      int[] fragmentEnds = new int[consensusBuilder.contigsFragmented.length];
      for(byte[] contig: consensusBuilder.contigsFragmented[0]){
        fragmentEnds[0] += contig.length;
      }
      for(int i = 1; i < consensusBuilder.contigsFragmented.length; i++){
        fragmentEnds[i] = fragmentEnds[i - 1];
        for(byte[] contig: consensusBuilder.contigsFragmented[i]){
          fragmentEnds[i] += contig.length;
        }
      }
      int fragmentID = 0;
      int contigID = 0;
      int predPos = 0;
      int assemblyPos = 0;
      for(RepeatInterval interval : repeatIntervals){
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
      for(RepeatInterval interval : repeatIntervals){
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
        if(interval.start - predPos > 0){
          contigs.add(Arrays.copyOfRange(finalConsensus, predPos, interval.start));
        }
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

  private int[] recomputeAssemblyConsensus(ConsensusBuilderWithReassembling consensusBuilder) throws InterruptedException{
//    Reference genome = new Reference(consensusBuilder.finalConsensus);
    byte[] consensus = consensusBuilder.getConsensusSimple(false, 0, consensusBuilder.finalConsensus.getBytes());
    consensusBuilder.finalConsensus = new String(consensus);
    int[] res = new int[consensus.length];
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
      ArrayList<String> contigs = new ArrayList<>();
      TByteArrayList consensusNew = new TByteArrayList();
      boolean isGap = true;
      int prev = 0;
      for(int i = 0, j = 0, k = 0; i < consensus.length; i++){
        if(consensus[i] == GAP){
          if(!isGap){
            // gap is started
            isGap = true;
            contigs.add(consensusBuilder.finalConsensus.substring(prev, i));
          }
          res[i] = -100000;
        }else{
          if(isGap){
            isGap = false;
            prev = i;
          }else if(i == consensusBuilder.contigEnds[j]){
            contigs.add(consensusBuilder.finalConsensus.substring(prev, i));
            prev = i;
          }
          consensusNew.add(consensus[i]);
          res[i] = k;
          k++;
        }
        if(i == consensusBuilder.contigEnds[j]){
          j++;
        }
      }
      contigs.add(consensusBuilder.finalConsensus.substring(prev));
      consensusBuilder.contigs = contigs.toArray(new String[contigs.size()]);
      consensusBuilder.contigEnds = new int[contigs.size()];
      for(int i = 0, j = 0; i < consensusBuilder.contigs.length; i++){
        j += consensusBuilder.contigs[i].length();
        consensusBuilder.contigEnds[i] = j;
      }
      consensusBuilder.finalConsensus = consensusNew.toString();
    }
    return res;
  }

  public void assemble() throws IOException, InterruptedException, InvocationTargetException, InstantiationException, IllegalAccessException, JDOMException{
    long time = System.currentTimeMillis();
    ArrayList<PairedRead> reads = dataReader.readFilesWithReads();
    ArrayList<Reference> selectedRefs = refFinder.selectReferences(reads);
    if(selectedRefs.size() == 0){
      Reference genome = new Reference(document);
      genome.reads = reads;
      genome.problemReads = new ArrayList<>();
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
      ProblemInterval[][] problemIntervals = new ProblemInterval[assemblies.length][];
      ConsensusBuilderWithReassembling[] cBuilders = new ConsensusBuilderWithReassembling[assemblies.length];
      int[] readNums = new int[assemblies.length];
      int[][] contigEnds = new int[assemblies.length][];
      for(int i = 0; i < selRefs.length; i++){
        ConsensusBuilderWithReassembling consensusBuilderGC = buildUsingRef(document, selRefs[i]);
        assemblies[i] = consensusBuilderGC.finalConsensus;
        mappedData[i] = consensusBuilderGC.mappedData;
        coordConvert[i] = consensusBuilderGC.genomeToAssembly;
        problemIntervals[i] = consensusBuilderGC.problemIntervals;
        cBuilders[i] = consensusBuilderGC;
        readNums[i] = consensusBuilderGC.reads.size();
        contigEnds[i] = consensusBuilderGC.contigEnds;
        Reference ref = selRefs[i];
        consensusBuilderGC.remapReadsToAssembly();
        consensusBuilderGC.printAssembly(ref.name + "_assembly_no_postproc.fasta");
      }
      ArrayList<RepeatInterval>[] intervals =
          postprocessor.findRegionsToCut(assemblies, contigEnds, selRefs, readNums);
      MappedData[] mappedDataAfterCuts = resortReads(intervals, mappedData);
      for(int i = 0; i < assemblies.length; i++){
        ConsensusBuilderWithReassembling cBuilder = cBuilders[i];
        reformatAssembly(cBuilder, intervals[i]);
        cBuilder.mappedData = mappedDataAfterCuts[i];
        for(MappedRead mappedRead : cBuilder.mappedData.mappedReads){
          adjustCoord(mappedRead, cBuilder.genomeToAssembly);
        }
        Reference ref = selRefs[i];
        cBuilder.mapper.mapReads(cBuilder.mappedData, new Reference(cBuilder));
        int[] coordNew = recomputeAssemblyConsensus(cBuilder);
        for(MappedRead mappedRead : cBuilder.mappedData.mappedReads){
          adjustCoord(mappedRead, coordNew);
        }
        cBuilder.printAssembly(ref.name + "_assembly.fasta");
        BamPrinter bamPrinter = new BamPrinter();
        Reference genome = new Reference(cBuilder);
        String bamPath = outPath + ref.name + "_mapped_reads.bam";
        bamPrinter.printBAM(cBuilder.mappedData, genome, bamPath);
        if(refFinder.debug && ref.problemReads.size() != 0){
          MappedData garbage = cBuilder.mapper.mapReads(ref.problemReads, genome);
          bamPath = outPath + ref.name + "_garbage.bam";
          bamPrinter.printBAM(garbage, genome, bamPath);
        }
      }
    }else{
      for(Reference ref : selectedRefs){
        buildUsingRefNoPostproc(document, ref);
      }
    }
    logger.println("Total time: " + (System.currentTimeMillis() - time)/1000);
  }

  private static void printUsage(){
    System.out.println("Usage: java -cp ./VirGenA.jar RefBasedAssembler pathToConfigFile");
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
