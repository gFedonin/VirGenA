/**
 * Created with IntelliJ IDEA.
 * Date: 04.08.16
 */
public class RandomModelBuilder{

  private static String pathToReference;
  private static String outPath = "./k5rn10000.rm";
  private static int scoreK = 5;
  private static int modelK = 0;
  private static float coef = 1.25f;
  private static int minReadLen = 10;
  private static int maxReadLen = 300;
  private static int step = 10;
  private static int readNum = 10000;
  private static int threadNum = 1;

  private static void printUsage(){
    System.out.println("Usage: java -cp ./VirGenA.jar RandomModelBuilder [Options]");
    System.out.println("Options:");
    System.out.println("[-ref|-refMSA] path to reference or reference MSA in fasta format, " +
        "mandatory");
    System.out.println("-k k-mer length, default 5");
    System.out.println("-order order of Markov model for random read generation, default k-1");
    System.out.println("-hitRatio hit concordance ratio used for k-mer score computation, " +
        "default 1.25");
    System.out.println("-readNum number of random reads of each length to use in the model, " +
        "default 10000");
    System.out.println("-minReadLen minimal length of random reads to use in the model, " +
        "default 10");
    System.out.println("-maxReadLen maximal length of random reads to use in the model, " +
        "default 300");
    System.out.println("-step number to increase random reads length starting from minReadLen " +
        "and ending with maxReadLen, should be between 1 and " +
        "(maxReadLen - minReadLen), default 10");
    System.out.println("-threadNum number of threads to use for computation, default 1");
    System.out.println("-out path to output file, default ./k5rn10000.rm");

  }

  public static void main(String[] args){
    boolean useMSA = false;
    if(args.length == 0){
      System.out.println("This program generates random reads of given lengths and computes " +
          "k-mer scores for these reads with given k and reference sequence or multiple " +
          "sequence alignment of reference sequences. Read lengths vary from minReadLen to " +
          "maxReadLen with given step. Defaults minReadLen=10, maxReadLen=300, step=10 mean" +
          "that random reads with lengths 10, 20, 30,...,290 and 300 will be generated. " +
          "The smaller the step - the more precise model you get, but more time it takes to " +
          "compute it. The bigger readsNum - the more reads will be generated and more precise " +
          "model you get.");
      printUsage();
      return;
    }
    for(int i = 0; i < args.length; i++){
      String arg = args[i];
      switch(arg){
        case "-ref":
          pathToReference = args[++i];
          break;
        case "-k":
          scoreK = Integer.parseInt(args[++i]);
          break;
        case "-hitRatio":
          coef = Float.parseFloat(args[++i]);
          break;
        case "-out":
          outPath = args[++i];
          break;
        case "-refMSA":
          useMSA = true;
          pathToReference = args[++i];
          break;
        case "-threadNum":
          threadNum = Integer.parseInt(args[++i]);
          break;
        case "-readNum":
          readNum = Integer.parseInt(args[++i]);
          break;
        case "-minReadLen":
          minReadLen = Integer.parseInt(args[++i]);
          break;
        case "-maxReadLen":
          maxReadLen = Integer.parseInt(args[++i]);
          break;
        case "-step":
          step = Integer.parseInt(args[++i]);
          break;
        case "-order":
          modelK = Integer.parseInt(args[++i]);
          break;
        default:
          System.out.printf("Unknown parameter name %s\n", arg);
          printUsage();
          return;
      }
    }
    if(modelK == 0){
      modelK = scoreK - 1;
    }
    if(step < 1 || step > maxReadLen - minReadLen){
      System.out.printf("Wrong step value: %d. Should be between 1 and " +
          "(maxReadLen - minReadLen)\n", step);
    }
    if(pathToReference == null){
      System.out.println("Path to reference or reference MSA must be provided!");
      printUsage();
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

}
