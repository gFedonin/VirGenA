import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;
import org.jdom2.Document;
import org.jdom2.Element;
import org.jdom2.input.SAXBuilder;

/**
 * Created with IntelliJ IDEA.
 * Date: 04.08.16
 */
public class RandomModelBuilder{

  static void addParameters(ArgumentParser parser){
    parser.description("This program " +
        "generates random reads of given lengths and computes k-mer scores for these reads with given k and reference " +
        "sequence or multiple sequence alignment of reference sequences. Read lengths vary from min_read_len to " +
        "max_read_len with given step. The smaller the step - the more precise model you get, but more time it takes to " +
        "compute it. The bigger read_num - the more reads will be generated and more precise model you get. The scores " +
        "obtained from random reads are used further to compute the threshold based on given p-value. For given p-value " +
        "you need at least 1/p-value reads. You should provide the path to config file, in which the parameter are " +
        "stored in <RandomModel> section. You can also specify some of the command line options explained below. Not " +
        "specified options will be taken from the config file.");
    parser.addArgument("-c").dest("config_file").help("Path to config file").required(true);
    parser.addArgument("-ref").dest("path_to_reference").help("Path to the reference or the reference MSA " +
        "in FASTA format.");
    parser.addArgument("--msa").action(Arguments.storeTrue()).help("Build MSA model, not single reference.");
    parser.addArgument("-k").help("Kmer length to use in Kmer score.");
    parser.addArgument("-order").dest("markov_model_order").help("The model order for random read " +
        "generation.");
    parser.addArgument("-indel").dest("indel_tolerance_threshold").help("The threshold ratio of distances " +
        "between two consequent k-mers in the read sequence to the corresponding distance in reference sequence. If the " +
        "computed distance is in [1/indel_tolerance_threshold, indel_tolerance_threshold] these k-mers are considered " +
        "to be concordant. Should be > 1.");
    parser.addArgument("-t").dest("thread_num").help("Number of threads to use.");
    parser.addArgument("-n").dest("read_num").help("Number of random reads to generate.");
    parser.addArgument("-minRL").dest("min_read_len").help("Minimal length of random reads to use in the model");
    parser.addArgument("-maxRL").dest("max_read_len").help("Maximal length of random reads to use in the model");
    parser.addArgument("-s").dest("step").help("number to increase random reads length starting from " +
        "min_read_len and ending with max_read_len, should be between 1 and (max_read_len - min_read_len)");
    parser.addArgument("-o").dest("out_path").help("path to output file.");
  }

  static void run(ArgumentParser parser, Namespace parsedArgs){
    boolean useMSA;
    SAXBuilder jdomBuilder = new SAXBuilder();
    Document jdomDocument;
    try{
      jdomDocument = jdomBuilder.build(parsedArgs.getString("config_file"));
    }catch(Exception e){
      e.printStackTrace();
      return;
    }
    Element root = jdomDocument.getRootElement();
    useMSA = parsedArgs.getBoolean("msa");
    String pathToReference = parsedArgs.getString("path_to_reference");
    String scoreKStr = parsedArgs.getString("k");
    String order = parsedArgs.getString("markov_model_order");
    String hitRatio = parsedArgs.getString("indel_tolerance_threshold");
    String tNum = parsedArgs.getString("thread_num");
    String rNum = parsedArgs.getString("read_num");
    String minRLen = parsedArgs.getString("min_read_len");
    String maxRLen = parsedArgs.getString("max_read_len");
    String s = parsedArgs.getString("step");
    String outPath = parsedArgs.getString("out_path");
    Element elem;
    int scoreK;
    if(useMSA){
      elem = root.getChild("ReferenceSelector").getChild("MapperToMSA").getChild("RandomModelParameters");
      if(pathToReference == null){
        pathToReference = root.getChild("ReferenceSelector").getChildText("ReferenceMSA");
      }
      if(scoreKStr == null){
        scoreK = Integer.parseInt(root.getChild("ReferenceSelector").getChild("MapperToMSA").getChildText("K"));
      }else{
        scoreK = Integer.parseInt(scoreKStr);
      }
      if(outPath == null){
        outPath = root.getChild("ReferenceSelector").getChild("MapperToMSA").getChildText("RandomModelPath");
      }
    }else{
      elem = root.getChild("Mapper").getChild("RandomModelParameters");
      if(pathToReference == null){
        pathToReference = root.getChildText("Reference");
      }
      if(scoreKStr == null){
        scoreK = Integer.parseInt(root.getChild("Mapper").getChildText("K"));
      }else{
        scoreK = Integer.parseInt(scoreKStr);
      }
      if(outPath == null){
        outPath = root.getChild("Mapper").getChildText("RandomModelPath");
      }
    }
    int modelK;
    if(order == null){
      modelK = Integer.parseInt(elem.getChildText("Order"));
    }else{
      modelK = Integer.parseInt(order);
    }
    float coef;
    if(hitRatio == null){
      coef = Float.parseFloat(elem.getChildText("IndelToleranceThreshold"));
    }else{
      coef = Float.parseFloat(hitRatio);
    }
    int threadNum;
    if(tNum == null){
      threadNum = Integer.parseInt(root.getChildText("ThreadNumber"));
    }else{
      threadNum = Integer.parseInt(tNum);
    }
    if(threadNum == -1){
      threadNum = Runtime.getRuntime().availableProcessors();
    }
    int minReadLen;
    if(minRLen == null){
      minReadLen = Integer.parseInt(elem.getChildText("MinReadLen"));
    }else{
      minReadLen = Integer.parseInt(minRLen);
    }
    int maxReadLen;
    if(maxRLen == null){
      maxReadLen = Integer.parseInt(elem.getChildText("MaxReadLen"));
    }else{
      maxReadLen = Integer.parseInt(maxRLen);
    }
    int readNum;
    if(rNum == null){
      readNum = Integer.parseInt(elem.getChildText("ReadNum"));
    }else{
      readNum = Integer.parseInt(rNum);
    }
    int step;
    if(s == null){
      step = Integer.parseInt(elem.getChildText("Step"));
    }else{
      step = Integer.parseInt(s);
    }
    if(step < 1 || step > maxReadLen - minReadLen){
      System.out.printf("Wrong step value: %d. Should be between 1 and " +
          "(maxReadLen - minReadLen)\n", step);
    }
    if(pathToReference.isEmpty()){
      System.out.println("Path to reference or reference MSA must be provided!");
      parser.printUsage();
      return;
    }
    try{
      if(useMSA){
        RandomModelMSA model = new RandomModelMSA(pathToReference, scoreK, modelK, coef, threadNum);
        model.printModel(readNum, minReadLen, maxReadLen, step, outPath);
      }else{
        RandomModel model = new RandomModel(pathToReference, scoreK, modelK, coef, threadNum);
        model.printModel(readNum, minReadLen, maxReadLen, step, outPath);
      }
    }catch(Exception e){
      e.printStackTrace();
    }
  }

  public static void main(String[] args){
    ArgumentParser parser = ArgumentParsers.newFor("RandomModelBuilder").build();
    addParameters(parser);
    if(args.length == 0){
      parser.printUsage();
      return;
    }
    Namespace parsedArgs;
    try{
      parsedArgs = parser.parseArgs(args);
      run(parser, parsedArgs);
    }catch(ArgumentParserException e){
      parser.handleError(e);
    }
  }

}
