/**
 * Created by Gennady
 */
class Hit extends Constants implements Comparable{
  int genomePos;
  int readPos;

  Hit(){

  }

  @Override
  public int compareTo(Object o){
    Hit h = (Hit) o;
    if(genomePos == h.genomePos){
      return readPos - h.readPos;
    }
    return genomePos - h.genomePos;
  }
}
