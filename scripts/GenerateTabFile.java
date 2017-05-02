import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.jdom2.Document;
import org.jdom2.input.SAXBuilder;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by Gennady on 02.05.2017.
 */
public class GenerateTabFile{

  private static String inPath;
  private static String outPath = "./";
  private static KMerCounter counter;
  private static Aligner aligner;

  private static int minReadLen = 50;
  private static byte[] complement = new byte[127];
  static {
    complement[65] = 84;
    complement[84] = 65;
    complement[71] = 67;
    complement[67] = 71;
    complement[45] = 45;
    complement[78] = 78;
  }

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
    MappedRead[] bestPositions;
    if(genome.isFragmented){
      bestPositions = counter.getNBestRegions(read, genome.index, genome.contigEnds);
    }else{
      bestPositions = counter.getNBestRegions(read, genome.index, 0, genome.length);
    }
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

  private static byte[] getComplement(byte[] seq) {
    int len = seq.length;
    byte[] res = new byte[len];

    for(int i = 0; i < len; ++i) {
      res[i] = complement[seq[len - i - 1]];
    }

    return res;
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
                writer.write(Integer.toString(mRead.start + 1) + "\t");
                writer.write(Integer.toString(mRead.end + 1) + "\t");
                writer.write(contigName1 + "\t");
                writer.write(Integer.toString(records[1].getAlignmentEnd()) + "\t");
                writer.write(Integer.toString(records[1].getAlignmentStart()) + "\n");
              }else{
                writer.write(contigName1 + "\t");
                writer.write(Integer.toString(records[1].getAlignmentStart()) + "\t");
                writer.write(Integer.toString(records[1].getAlignmentEnd()) + "\t");
                writer.write(contigName0 + "\t");
                writer.write(Integer.toString(mRead.end + 1) + "\t");
                writer.write(Integer.toString(mRead.start + 1) + "\n");
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
                writer.write(Integer.toString(mRead.start + 1) + "\t");
                writer.write(Integer.toString(mRead.end + 1) + "\t");
                writer.write(contigName0 + "\t");
                writer.write(Integer.toString(records[0].getAlignmentEnd()) + "\t");
                writer.write(Integer.toString(records[0].getAlignmentStart()) + "\n");
              }else{
                writer.write(contigName0 + "\t");
                writer.write(Integer.toString(records[0].getAlignmentStart()) + "\t");
                writer.write(Integer.toString(records[0].getAlignmentEnd()) + "\t");
                writer.write(contigName1 + "\t");
                writer.write(Integer.toString(mRead.end + 1) + "\t");
                writer.write(Integer.toString(mRead.start + 1) + "\n");
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
              writer.write(Integer.toString(records[1].getAlignmentStart()) + "\t");
              writer.write(Integer.toString(records[1].getAlignmentEnd()) + "\t");
              writer.write(contigName0 + "\t");
              writer.write(Integer.toString(records[0].getAlignmentEnd()) + "\t");
              writer.write(Integer.toString(records[0].getAlignmentStart()) + "\n");
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
              writer.write(Integer.toString(records[0].getAlignmentStart()) + "\t");
              writer.write(Integer.toString(records[0].getAlignmentEnd()) + "\t");
              writer.write(contigName1 + "\t");
              writer.write(Integer.toString(records[1].getAlignmentEnd()) + "\t");
              writer.write(Integer.toString(records[1].getAlignmentStart()) + "\n");
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
        writerLib.write(sampleName + " TAB " + tabPath + " " +
            Integer.toString(avg) + " 0.01 FR");
      }else{
        writerLib.write(sampleName + " TAB " + tabPath + " " +
            Integer.toString(avg) + " " + Float.toString((float) sigma/avg) + " FR");
      }
    }else{
      writerLib.write(sampleName + " TAB " + tabPath + " " +
          Integer.toString(avg) + " 0.99 FR");
    }
    writerLib.close();
  }

  private static HashMap<String, String> readFasta(String path) throws IOException{
    HashMap<String, String> res = new HashMap<>();
    BufferedReader reader = new BufferedReader(new FileReader(path));
    StringBuilder builder = new StringBuilder();
    String name = null;
    for(String line = reader.readLine(); line != null; line = reader.readLine()){
      if(line.startsWith(">")){
        if(name != null){
          res.put(name, builder.toString());
        }
        name = line.substring(1);
        builder = new StringBuilder();
      }else{
        builder.append(line);
      }
    }
    if(name != null){
      res.put(name, builder.toString());
    }
    return res;
  }

  private static void printUsage(){
    System.out.println("Usage: java -cp ./VirGenA.jar:. GenerateTabFile [Options]");
    System.out.println("Options:");
    System.out.println("-in path to input folder which contains VirGenA assembly with BAM " +
        "files, generated by VirGenA. This is a mandatory parameter!");
    System.out.println("-out path to output folder to which TAB and LIB files will be printed, default './'.");
    System.out.println("-c path to a valid VirGenA config file: please use the one you " +
        "used for assembly. " +
        "This script will use the same random model file and the other Mapper parameters to " +
        "align reads to the assembly. This is a mandatory parameter!");
    System.out.println("-rl minimal read length to use for scaffolding, " +
        "default 50");

  }

  private static File parseArgs(String[] args){
    if(args.length == 0 || args[0].equals("-h") || args[0].equals("--help")){
      System.out.println("This program generates TAB file and corresponding LIB file " +
          " needed to run SSPACE scaffolder using VirGenA assembly with BAM file.");
      printUsage();
      return null;
    }
    File inFolder = null;
    File outFolder;
    for(int i = 0; i < args.length; i+=2){
      switch(args[i]){
        case "-in":
          inPath = args[i + 1];
          inFolder = new File(inPath);
          if(!inFolder.exists()){
            System.out.println("No such folder: " + inPath);
            return null;
          }else{
            if(!inFolder.isDirectory()){
              System.out.println(inPath + " should be an input folder with VirGenA assembly.");
              return null;
            }
          }
          break;
        case "-out":
          outPath = args[i + 1];
          outFolder = new File(outPath);
          if(!outFolder.exists()){
            outFolder.mkdirs();
          }
          break;
        case "-c":
          String configPath = args[i + 1];
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
          break;
        case "-rl":
          minReadLen = Integer.parseInt(args[i + 1]);
          break;
      }
    }
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

  public static void main(String[] args){
    try{
      File inFolder = parseArgs(args);
      if(inFolder == null){
        return;
      }
      String sampleName = inFolder.getName();
      ArrayList<String> bamFiles = new ArrayList<>();
      HashMap<String, String> contigs = new HashMap<>();
      for(File file : inFolder.listFiles()){
        if(file.getName().endsWith("_assembly.fasta")){
          contigs.putAll(readFasta(file.getPath()));
        }else if(file.getName().endsWith(".bam")){
          bamFiles.add(file.getPath());
        }
      }
      if(contigs.size() == 1){
        System.out.println("There is only one contig: nothing to scaffold!");
        return;
      }
      printTAB(bamFiles, contigs, sampleName);
      BufferedWriter writer = new BufferedWriter(new FileWriter(outPath +
          sampleName + "_contigs.fasta"));
      for(HashMap.Entry<String, String> contig: contigs.entrySet()){
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

}