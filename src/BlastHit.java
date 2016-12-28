import java.util.ArrayList;

/**
 * Created by Gennady on 02.04.2016.
 */
public class BlastHit implements Comparable{

  public String hitID;
  public int queryLen;
  public int match;
  public ArrayList<HSP> hsps;
  public int bestHspLength;
  public float bestHspCoverage;
  public float bestHspNormScore;
  public HSP bestHsp;
  public float totalScore;

  public BlastHit(){
    hsps = new ArrayList<>();
  }

  public void computeStats(){
    if(hsps.size() == 0){
      return;
    }
    bestHsp = hsps.get(0);
    float mScore = bestHsp.score;
    for(HSP hsp: hsps){
      if(hsp.score > mScore){
        bestHsp = hsp;
        mScore = hsp.score;
      }
      totalScore += hsp.score;
    }
    bestHspLength = bestHsp.endQ - bestHsp.startQ;
    bestHspCoverage = (float)bestHspLength/queryLen;
    bestHspNormScore = bestHsp.score/match/queryLen;
  }

  @Override
  public int compareTo(Object o){
    BlastHit hit = (BlastHit) o;
//    if(hit.avIdentity == avIdentity){
//      return 0;
//    }
//    return (hit.avIdentity > avIdentity) ? 1 : -1;
//    if(hit.avScore == avScore){
//      return 0;
//    }
//    return (hit.avScore > avScore) ? 1 : -1;
    if(hit.totalScore == totalScore){
      return 0;
    }
    return (hit.totalScore > totalScore) ? 1 : -1;
  }

}
