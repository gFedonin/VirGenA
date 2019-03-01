import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.Namespace;
import org.jdom2.Document;
import org.jdom2.input.SAXBuilder;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by Gennady on 02.05.2017.
 */
public class GenerateTabFile extends Constants{

  private static String inPath;
  private static String outPath = "./";
  private static KMerCounter counter;
  private static Aligner aligner;

  private static int minReadLen = 50;

  private static Reference makeReference(HashMap<String, String> contigs){
    String[] contigNames = new String[contigs.size()];
    int[] ends = new int[contigNames.length];
    int i = 0;
    StringBuilder builder = new StringBuilder();
    for(Map.Entry<String, String> contig: contigs.entrySet()){
      contigNames[i] = contig.getKey();
      String seq = contig.getValue();
      builder.append(seq);
      ends[i] = builder.length();
      i++;
    }
    int threadNum = Runtime.getRuntime().availableProcessors();
    return new Reference("", builder.toString(), counter.K, ends, ends,
        contigNames, true, threadNum);
  }

  private static MappedRead[] mapReadFast(byte[] read, Reference genome){
    MappedRead[] bestPositions = counter.getNBestRegions(read, genome.index, genome.contigEnds);

    if(bestPositions == null){
      return new MappedRead[0];
    }
    for(MappedRead mappedRead : bestPositions){
      mappedRead.seq = read;
      mappedRead.aln = aligner.align(read, Arrays.copyOfRange(genome.seqB,
          mappedRead.start, mappedRead.end));
    }
    return bestPositions;
  }

  private static String fixContigName(String contigName, HashMap<String, String> contigs){
    if(!contigs.containsKey(contigName)){
      boolean found = false;
      for(String cName : contigs.keySet()){
        if(cName.contains(contigName)){
          contigName = cName;
          found = true;
          break;
        }
      }
      if(!found){
        return null;
      }
    }
    return contigName;
  }

  private static String getContigName(MappedRead mRead, Reference ref){
    for(int i = 0; i < ref.fragmentEnds.length; i++){
      if(mRead.start + mRead.aln.end2 <= ref.fragmentEnds[i]){
        if(i > 0){
          mRead.end = mRead.start + mRead.aln.end2 - ref.fragmentEnds[i - 1];
          mRead.start += mRead.aln.start2 - ref.fragmentEnds[i - 1];
        }else{
          mRead.end = mRead.start + mRead.aln.end2;
          mRead.start += mRead.aln.start2;
        }
        mRead.aln = null;
        return ref.fragmentNames[i];
      }
    }
    return null;
  }

  private static void printTAB(ArrayList<String> bamFiles, HashMap<String, String> contigs,
                               String sampleName) throws IOException{
    int sumInsertSize = 0;
    long sumInsSizeSqr = 0;
    int pairNum = 0;
    Reference ref = makeReference(contigs);
    File tabFile = new File(outPath + sampleName + ".tab");
    BufferedWriter writer = new BufferedWriter(new FileWriter(tabFile));
    for(String bamPath: bamFiles){
      SamReaderFactory factory = SamReaderFactory.makeDefault();
      factory.validationStringency(ValidationStringency.SILENT);
      SamReader sReader = factory.open(new File(bamPath));
      HashMap<String, SAMRecord[]> readToMappings = new HashMap<>();
      for(SAMRecord record : sReader){
        if(record.getReadLength() < minReadLen){
          continue;
        }
        SAMRecord[] records = readToMappings.get(record.getReadName());
        if(records == null){
          records = new SAMRecord[2];
          readToMappings.put(record.getReadName(), records);
        }
        if(record.getFirstOfPairFlag()){
          records[0] = record;
        }else{
          records[1] = record;
        }
      }
      for(SAMRecord[] records: readToMappings.values()){
        if(records[0] == null || records[1] == null){
          continue;
        }
        if(records[0].getReadUnmappedFlag()){
          if(!records[1].getReadUnmappedFlag()){
            String contigName1 = records[1].getReferenceName();
            contigName1 = fixContigName(contigName1, contigs);
            MappedRead[] mReads;
            if(records[1].getReadNegativeStrandFlag()){
              mReads = mapReadFast(records[0].getReadBases(), ref);
            }else{
              mReads = mapReadFast(getComplement(records[0].getReadBases()), ref);
            }
            boolean[] keep = new boolean[mReads.length];
            int maxScore = 0;
            for(MappedRead mRead : mReads){
              if(mRead.aln.score > maxScore){
                maxScore = mRead.aln.score;
              }
            }
            for(int i = 0; i < mReads.length; i++){
              if(mReads[i].aln.score == maxScore){
                keep[i] = true;
              }
            }
            for(int i = 0; i < mReads.length; i++){
              MappedRead mRead = mReads[i];
              if(!keep[i]){
                continue;
              }
              String contigName0 = getContigName(mRead, ref);
              if(records[1].getReadNegativeStrandFlag()){
                writer.write(contigName0 + "\t");
                writer.write((mRead.start + 1) + "\t");
                writer.write((mRead.end + 1) + "\t");
                writer.write(contigName1 + "\t");
                writer.write(records[1].getAlignmentEnd() + "\t");
                writer.write(records[1].getAlignmentStart() + "\n");
              }else{
                writer.write(contigName1 + "\t");
                writer.write(records[1].getAlignmentStart() + "\t");
                writer.write(records[1].getAlignmentEnd() + "\t");
                writer.write(contigName0 + "\t");
                writer.write((mRead.end + 1) + "\t");
                writer.write((mRead.start + 1) + "\n");
              }
            }
          }
        }else{
          String contigName0 = records[0].getReferenceName();
          contigName0 = fixContigName(contigName0, contigs);
          if(records[1].getReadUnmappedFlag()){
            MappedRead[] mReads;
            if(records[0].getReadNegativeStrandFlag()){
              mReads = mapReadFast(records[1].getReadBases(), ref);
            }else{
              mReads = mapReadFast(getComplement(records[1].getReadBases()), ref);
            }
            boolean[] keep = new boolean[mReads.length];
            int maxScore = 0;
            for(MappedRead mRead : mReads){
              if(mRead.aln.score > maxScore){
                maxScore = mRead.aln.score;
              }
            }
            for(int i = 0; i < mReads.length; i++){
              if(mReads[i].aln.score == maxScore){
                keep[i] = true;
              }
            }
            for(int i = 0; i < mReads.length; i++){
              MappedRead mRead = mReads[i];
              if(!keep[i]){
                continue;
              }
              String contigName1 = getContigName(mRead, ref);
              if(records[0].getReadNegativeStrandFlag()){
                writer.write(contigName1 + "\t");
                writer.write((mRead.start + 1) + "\t");
                writer.write((mRead.end + 1) + "\t");
                writer.write(contigName0 + "\t");
                writer.write((records[0].getAlignmentEnd()) + "\t");
                writer.write((records[0].getAlignmentStart()) + "\n");
              }else{
                writer.write(contigName0 + "\t");
                writer.write((records[0].getAlignmentStart()) + "\t");
                writer.write((records[0].getAlignmentEnd()) + "\t");
                writer.write(contigName1 + "\t");
                writer.write((mRead.end + 1) + "\t");
                writer.write((mRead.start + 1) + "\n");
              }
            }
          }else{
            String contigName1 = records[1].getReferenceName();
            contigName1 = fixContigName(contigName1, contigs);
            if(records[0].getReadNegativeStrandFlag()){
              if(records[1].getReadNegativeStrandFlag()){
                continue;
              }
              if(contigName0.equals(contigName1)){
                int size = records[0].getAlignmentEnd() - records[1].getAlignmentStart();
                sumInsertSize += size;
                sumInsSizeSqr += size*size;
                pairNum++;
              }
              writer.write(contigName1 + "\t");
              writer.write(records[1].getAlignmentStart() + "\t");
              writer.write(records[1].getAlignmentEnd() + "\t");
              writer.write(contigName0 + "\t");
              writer.write(records[0].getAlignmentEnd() + "\t");
              writer.write(records[0].getAlignmentStart() + "\n");
            }else{
              if(!records[1].getReadNegativeStrandFlag()){
                continue;
              }
              if(contigName0.equals(contigName1)){
                int size = records[1].getAlignmentEnd() - records[0].getAlignmentStart();
                sumInsertSize += size;
                sumInsSizeSqr += size*size;
                pairNum++;
              }
              writer.write(contigName0 + "\t");
              writer.write(records[0].getAlignmentStart() + "\t");
              writer.write(records[0].getAlignmentEnd() + "\t");
              writer.write(contigName1 + "\t");
              writer.write(records[1].getAlignmentEnd() + "\t");
              writer.write(records[1].getAlignmentStart() + "\n");
            }
          }
        }
      }
      sReader.close();
    }
    writer.close();
    int avg = sumInsertSize/pairNum;
    long sigma = Math.round(Math.sqrt(sumInsSizeSqr/pairNum - avg*avg));
    BufferedWriter writerLib = new BufferedWriter(new FileWriter(outPath + sampleName + ".lib"));
    String tabPath = tabFile.getCanonicalPath();
    if(sigma < avg){
      if(sigma == 0){
        writerLib.write(sampleName + " TAB " + tabPath + " " + avg + " 0.01 FR");
      }else{
        writerLib.write(sampleName + " TAB " + tabPath + " " + avg + " " + (float) sigma/avg + " FR");
      }
    }else{
      writerLib.write(sampleName + " TAB " + tabPath + " " + avg + " 0.99 FR");
    }
    writerLib.close();
  }

  private static File init(Namespace parsedArgs){
    inPath = parsedArgs.getString("in_path");
    File inFolder = new File(inPath);
    if(!inFolder.exists()){
      System.out.println("No such folder: " + inPath);
      return null;
    }else{
      if(!inFolder.isDirectory()){
        System.out.println(inPath + " should be an input folder with VirGenA assembly.");
        return null;
      }
    }
    outPath = parsedArgs.getString("out_path");
    File outFolder = new File(outPath);
    if(!outFolder.exists()){
      outFolder.mkdirs();
    }
    String configPath = parsedArgs.getString("config_file");
    File configFile = new File(configPath);
    if(!configFile.exists()){
      System.out.println("No such file: " + configPath);
      return null;
    }
    SAXBuilder jdomBuilder = new SAXBuilder();
    Document jdomDocument;
    try{
      jdomDocument = jdomBuilder.build(configPath);
    }catch(Exception e){
      System.out.println("Wrong config file format! Please use valid VirGenA config file.");
      return null;
    }
    counter = KMerCounter.getInstance(jdomDocument);
    aligner = new SmithWatermanGotoh(jdomDocument);
    minReadLen = Integer.parseInt(parsedArgs.getString("min_read_len"));

    if(inPath == null){
      System.out.println("Input path is missing! Please provide it with -in option!");
      return null;
    }
    if(counter == null){
      System.out.println("Config file is missing! Please provide it with -c option!");
      return null;
    }
    return inFolder;
  }

  static void run(Namespace parsedArgs){
    File inFolder = init(parsedArgs);
    String sampleName = inFolder.getName();
    ArrayList<String> bamFiles = new ArrayList<>();
    HashMap<String, String> contigs = new HashMap<>();
    for(File file : inFolder.listFiles()){
      if(file.getName().endsWith("_assembly.fasta")){
        try{
          contigs.putAll(FastaReader.read(file.getPath()));
        }catch(Exception e){
          e.printStackTrace();
        }
      }else if(file.getName().endsWith(".bam")){
        bamFiles.add(file.getPath());
      }
    }
    if(contigs.size() == 1){
      System.out.println("There is only one contig: nothing to scaffold!");
      return;
    }
    try{
      printTAB(bamFiles, contigs, sampleName);
      BufferedWriter writer = new BufferedWriter(new FileWriter(outPath +
          sampleName + "_contigs.fasta"));
      for(HashMap.Entry<String, String> contig : contigs.entrySet()){
        if(contig.getValue().length() > 0){
          writer.write(">" + contig.getKey() + "\n");
          writer.write(contig.getValue() + "\n");
        }
      }
      writer.close();
    }catch(Exception e){
      e.printStackTrace();
    }
  }

  static void addParameters(ArgumentParser parser){
    parser.description("This program generates TAB file and corresponding LIB file " +
        "needed to run SSPACE scaffolder using VirGenA assembly with BAM file.");
    parser.addArgument("-c").dest("config_file").help("path to a valid VirGenA config file: please use the " +
        "one you used for assembly.").required(true);
    parser.addArgument("-in").dest("in_path").help("path to a valid VirGenA config file: please use the " +
        "one you used for assembly.").required(true);
    parser.addArgument("-out").dest("out_path").help("path to output folder to which TAB and LIB files will " +
        "be printed, default './'.").required(false);
    parser.addArgument("-rl").dest("min_read_len").help("minimal read length to use for scaffolding, " +
        "default 50").required(false);
  }

  public static void main(String[] args){
    try{
      ArgumentParser parser = ArgumentParsers.newFor("GenerateTabFile").build();
      addParameters(parser);
      if(args.length == 0){
        parser.printUsage();
        return;
      }
      run(parser.parseArgs(args));
    }catch(Exception e){
      e.printStackTrace();
    }
  }

}