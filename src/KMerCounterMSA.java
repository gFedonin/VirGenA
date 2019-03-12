import org.jdom2.Document;
import org.jdom2.Element;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * Date: 03.01.15
 */
class KMerCounterMSA extends KMerCounterBase{

  private static KMerCounterMSA instance;

  private KMerCounterMSA(Document document){
    try {
      Element element = document.getRootElement().getChild("ReferenceSelector").getChild("MapperToMSA");
      pValue = Float.parseFloat(element.getChildText("pValue"));
      DataReader dataReader = DataReader.getInstance(document);
//      randomReadsCountsPath = element.getChildText("RandomModelPath");
//      if(randomReadsCountsPath.isEmpty()){
        K = Integer.parseInt(element.getChildText("K"));
        highCoef = Float.parseFloat(element.getChildText("IndelToleranceThreshold"));
        lowCoef = 1/highCoef;
        element = element.getChild("RandomModelParameters");
//        String outPath = document.getRootElement().getChildText("OutPath");
        Logger logger = Logger.getInstance(document);
//        logger.println("No path to random read model found in the config file for <MapperToMSA>. Creating new model " +
//            "with the parameters given in the config file in <RandomModelParameters> and K = " + scoreK);
      logger.println("Creating random reads model from MSA with the parameters given in the config file in <RandomModelParameters>");
        RandomModelMSA model = new RandomModelMSA(document);
//        int minReadLen = Integer.parseInt(element.getChildText("MinReadLen"));
//        int maxReadLen = Integer.parseInt(element.getChildText("MaxReadLen"));
        int step = Integer.parseInt(element.getChildText("Step"));
        int readNum = Integer.parseInt(element.getChildText("ReadNum"));
//        randomReadsCountsPath = outPath + "random_model_msa.rm";
//        model.printModel(readNum, minReadLen, maxReadLen, step, randomReadsCountsPath);
//        logger.println("Random model is created and saved in " + randomReadsCountsPath);
//        logger.println("You can reuse it in future by providing the path to random_model.rm in <RandomModelPath> in " +
//            "config file in <MapperToMSA> section. In case then such path is provided, the model parameters in " +
//            "<RandomModelParameters> section and the K in <MapperToMSA> section are ignored.");
      model.genModel(readNum, dataReader.minReadLen, dataReader.maxReadLen, step);
//      }
      computeCuttofs(dataReader.minReadLen, dataReader.maxReadLen, step, model.counts);
//      readRandomModel();
    }catch (Exception e){
      e.printStackTrace();
    }
  }

  private KMerCounterMSA(int K, float coef){
    this.K = K;
    this.highCoef = coef;
    this.lowCoef = 1/coef;
  }

  static KMerCounterMSA getInstance(Document document){
    if(instance == null){
      instance = new KMerCounterMSA(document);
    }
    return instance;
  }

  static KMerCounterMSA getInstance(int K, float coef){
    if(instance == null){
      instance = new KMerCounterMSA(K, coef);
    }
    return instance;
  }


  MappedRead[] getNBestRegions(byte[] read, HashMap<String, int[]> index, int from, int to, int minCount){
    String readStr = new String(read);
    LongHit[] hits = getHits(readStr, index);
    if(hits.length == 0){
      return null;
    }
    int readLen = read.length;
    int totalHit = hits.length;
    boolean[] isConcordant = concordanceArray(hits);

    Location location = getInitState(from, to, readLen, hits, isConcordant);
    ArrayList<MappedRead> bestPositions = new ArrayList<>();
    if(location.count > minCount){
      bestPositions.add(new MappedRead(location.start, location.end, location.count));
    }
    for(int i = 1; i < totalHit; i++){
      computeCount(from, to, readLen, hits, i, location, isConcordant);
      if(location.count > minCount){
        bestPositions.add(new MappedRead(location.start, location.end, location.count));
      }
    }
    if(bestPositions.size() == 0){
      return null;
    }

    return pruneWindows(bestPositions);
  }

}
