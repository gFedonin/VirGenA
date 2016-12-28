/**
 * Created with IntelliJ IDEA.
 * Date: 20.01.14
 */
class LongHit extends Hit{

  public int len;

  LongHit(int gPos, int rPos){
    genomePos = gPos;
    readPos = rPos;
  }

  public boolean equals(Object o){
    LongHit hit = (LongHit) o;
    return (genomePos == hit.genomePos) && (readPos == hit.readPos) && (len == hit.len);
  }

}
