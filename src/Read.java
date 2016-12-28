/**
 * Created by Gennady
 */
class Read implements Comparable{

  String name;
  byte[] seq;
  byte[] q;
  char n;

  Read(String name, byte[] seq, char n){
    this.name = name;
    this.seq = seq;
    this.n = n;
  }

  @Override
  public int compareTo(Object o){
    Read r = (Read) o;
    return name.compareTo(r.name);
  }

}
