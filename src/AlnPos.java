/**
 * Created with IntelliJ IDEA.
 * Date: 03.01.15
 */
class AlnPos implements Comparable{
  short pos;
  short[] refIDs;

  @Override
  public int compareTo(Object o){
    AlnPos alnPos = (AlnPos) o;
    return pos - alnPos.pos;
  }
}
