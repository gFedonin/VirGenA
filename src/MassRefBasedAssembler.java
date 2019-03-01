import org.jdom2.Document;
import org.jdom2.input.SAXBuilder;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by Gennady on 21.02.2016.
 */
public class MassRefBasedAssembler{

  private static void printUsage(){
    System.out.println("Usage: java -cp ./ViRGA.jar MassRefBasedAssembler [Options]");
    System.out.println("Options:");
    System.out.println("-config path to config file");
    System.out.println("-in path to folder with pairedReads, default ./");
    System.out.println("-out path to output folder, default ./ref_asm/");
    System.out.println("-list path to the list of samples - a subset of files in the input " +
        "folder, one sample ID per file line, each sample ID should be the prefix of two file names " +
        "in the input folder. If passed, then -suffix option is ignored.");
    System.out.println("-suffix1,2 substrings containing suffixes of all files in input folder, " +
        "that should be processes. These files will be grouped by file name prefixes, i.e. by " +
        "substrings obtained by trimming given suffixes.");
    System.out.println("-overwrite if not set - for each sample checks if output dir exists. " +
        "If so, checks if assembly is done well. If so - do nothing, else - reassemble. " +
        "If set - reassemble anyway.");
  }

  public static void main(String[] args){
    try{
      long time = System.currentTimeMillis();
      if(args.length == 0){
        printUsage();
        return;
      }
      String pathToConfig = "config.xml";
      String pathToInFolder = "./";
      String pathToOutFolder = "./ref_asm/";
      String pathToList = null;
      String suffix1 = "_1.fastq.gz";
      String suffix2 = "_2.fastq.gz";
      boolean reassemble = false;
      for(int i = 0; i < args.length; i++){
        String arg = args[i];
        switch(arg){
          case "-config":
            pathToConfig = args[++i];
            break;
          case "-out":
            pathToOutFolder = args[++i];
            break;
          case "-in":
            pathToInFolder = args[++i];
            break;
          case "-list":
            pathToList = args[++i];
            break;
          case "-suffix1":
            suffix1 = args[++i];
            break;
          case "-suffix2":
            suffix2 = args[++i];
            break;
          case "-overwrite":
            reassemble = true;
            break;
          default:
            System.out.printf("Unknown parameter name %s\n", arg);
            printUsage();
            return;
        }
      }
      SAXBuilder jdomBuilder = new SAXBuilder();
      Document config = jdomBuilder.build(pathToConfig);
      File inFolder = new File(pathToInFolder);
      if(!inFolder.exists() || !inFolder.isDirectory()){
        System.err.printf("No such directory %s", pathToInFolder);
        return;
      }
      File outFolder = new File(pathToOutFolder);
      if(!outFolder.exists()){
        outFolder.mkdirs();
      }
      HashMap<String, ArrayList<String>> sampleIDsToFiles = new HashMap<>();
      if(pathToList != null){
        ArrayList<String> sampleIDs = new ArrayList<>();
        BufferedReader reader = new BufferedReader(new FileReader(pathToList));
        for(String line = reader.readLine(); line != null; line = reader.readLine()){
          sampleIDs.add(line);
        }
        for(File file: inFolder.listFiles()){
          String fileName = file.getName();
          for(String sampleID : sampleIDs){
            if(fileName.startsWith(sampleID)){
              ArrayList<String> fileNames = sampleIDsToFiles.get(sampleID);
              if(fileNames == null){
                fileNames = new ArrayList<>();
                sampleIDsToFiles.put(sampleID, fileNames);
              }
              fileNames.add(fileName);
            }
          }
        }
      }else{
        for(File file: inFolder.listFiles()){
          String fileName = file.getName();
          int i = fileName.indexOf(suffix1);
          if(i == -1){
            i = fileName.indexOf(suffix2);
          }
          if(i == -1){
            continue;
          }
          String sampleID = fileName.substring(0, i);
          ArrayList<String> fileNames = sampleIDsToFiles.get(sampleID);
          if(fileNames == null){
            fileNames = new ArrayList<>();
            sampleIDsToFiles.put(sampleID, fileNames);
          }
          fileNames.add(fileName);
        }
      }
      for(Map.Entry<String, ArrayList<String>> entry: sampleIDsToFiles.entrySet()){
        ArrayList<String> fileNames = entry.getValue();
        Collections.sort(fileNames);
        String sampleID = entry.getKey();
        if(fileNames.size() == 1){
          System.err.printf("Sample %s. Paired reads should be in two different files! " +
              "This ID will be ignored!\n", sampleID);
          continue;
        }
        if(fileNames.size() > 2){
          System.err.printf("Sample %s. Too many files with same sampleID! " +
              "This ID will be ignored!\n", sampleID);
          continue;
        }

        config.getRootElement().getChild("Data").getChild("pathToReads1").setText(pathToInFolder +
            fileNames.get(0));
        config.getRootElement().getChild("Data").getChild("pathToReads2").setText(pathToInFolder +
            fileNames.get(1));
        config.getRootElement().getChild("OutPath").setText(pathToOutFolder +
            sampleID + "/");
        boolean selectRefs = Boolean.parseBoolean(
            config.getRootElement().getChild("ReferenceSelector").getChildText("Enabled"));
        if(selectRefs){
          File file = new File(pathToOutFolder +
              sampleID + "/");
          if(reassemble){
            System.out.println("Working with " + sampleID);
            RefBasedAssembler assembler = new RefBasedAssembler(config);
            assembler.assemble();
          }else{
            if(file.exists()){
              File logFile = new File(pathToOutFolder +
                  sampleID + "/log.txt");
              boolean repeat = false;
              if(logFile.exists()){
                BufferedReader reader = new BufferedReader(new FileReader(logFile));
                int refNum = 0;
                for(String line = reader.readLine(); line != null; line = reader.readLine()){
                  if(line.startsWith("Selected references")){
                    line = reader.readLine();
                    while(!line.startsWith("Reference selection time")){
                      refNum++;
                      line = reader.readLine();
                      if(line == null){
                        repeat = true;
                        break;
                      }
                    }
                    break;
                  }
                }
                reader.close();
                if(!repeat){
                  int assemblyDone = 0;
                  int alignmentDone = 0;
                  for(File f : file.listFiles()){
                    if(f.getName().endsWith("_assembly.fasta")){
                      assemblyDone++;
                    }else if(f.getName().endsWith(".bam")){
                      alignmentDone++;
                    }
                  }
                  if(assemblyDone == 0 || assemblyDone < refNum || alignmentDone < refNum){
                    repeat = true;
                  }
                }
              }else{
                repeat = true;
              }
              if(repeat){
                System.out.println("Working with " + sampleID);
                RefBasedAssembler assembler = new RefBasedAssembler(config);
                try{
                  assembler.assemble();
                }catch(Exception e){
                  e.printStackTrace(System.err);
                  System.out.println("Problems with " + sampleID + ". Will keep working with the rest.");
                }
              }
            }else{
              System.out.println("Working with " + sampleID);
              RefBasedAssembler assembler = new RefBasedAssembler(config);
              try{
                assembler.assemble();
              }catch(Exception e){
                e.printStackTrace(System.err);
                System.out.println("Problems with " + sampleID + ". Will keep working with the rest.");
              }
            }
          }
        }else{
          File file = new File(pathToOutFolder +
              sampleID + "/");
          if(reassemble){
            System.out.println("Working with " + sampleID);
            ConsensusBuilderWithReassembling cBuilder = new ConsensusBuilderWithReassembling(config);
            try{
              cBuilder.assemble(config);
            }catch(Exception e){
              e.printStackTrace(System.err);
              System.out.println("Problems with " + sampleID + ". Will keep working with the rest.");
            }
          }else{
            if(file.exists()){
              File logFile = new File(pathToOutFolder +
                  sampleID + "/log.txt");
              if(logFile.exists()){
                int assemblyDone = 0;
                int alignmentDone = 0;
                for(File f : file.listFiles()){
                  if(f.getName().equals("assembly.fasta")){
                    assemblyDone++;
                  }else if(f.getName().endsWith(".bam")){
                    alignmentDone++;
                  }
                }
                if(assemblyDone != 1 || alignmentDone != 1){
                  System.out.println("Working with " + sampleID);
                  ConsensusBuilderWithReassembling cBuilder = new ConsensusBuilderWithReassembling(config);
                  try{
                    cBuilder.assemble(config);
                  }catch(Exception e){
                    e.printStackTrace(System.err);
                    System.out.println("Problems with " + sampleID + ". Will keep working with the rest.");
                  }
                }
              }else{
                System.out.println("Working with " + sampleID);
                ConsensusBuilderWithReassembling cBuilder = new ConsensusBuilderWithReassembling(config);
                try{
                  cBuilder.assemble(config);
                }catch(Exception e){
                  e.printStackTrace(System.err);
                  System.out.println("Problems with " + sampleID + ". Will keep working with the rest.");
                }
              }
            }else{
              System.out.println("Working with " + sampleID);
              ConsensusBuilderWithReassembling cBuilder = new ConsensusBuilderWithReassembling(config);
              try{
                cBuilder.assemble(config);
              }catch(Exception e){
                e.printStackTrace(System.err);
                System.out.println("Problems with " + sampleID + ". Will keep working with the rest.");
              }
            }
          }
        }
      }
      System.out.printf("Total time: %d, s", (System.currentTimeMillis() - time)/1000);
    }catch(Exception e){
      e.printStackTrace();
    }
  }

}
