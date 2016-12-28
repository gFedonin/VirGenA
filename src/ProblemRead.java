import java.util.Arrays;

/**
 * Created by Геннадий on 06.08.2014.
 */
class ProblemRead extends MappedRead{

  short gapID;
  EdgeAlignment aln;
  short lastLen;
  short countLeft;
  short lastLenLeft;
  short countRight;
  short lastLenRight;
  byte mappingType; // 0 left, 1 right, 2 both, -1 unmapped
  short cutoff;
  MappedRead mate;

  ProblemRead(String name, byte n, byte[] seq, short gapID, byte reverse){
    this.name = name;
    this.seq = seq;
    this.gapID = gapID;
    this.reverse = reverse;
    this.n = n;
    mappingType = -1;
  }


  void setMapping(MappedRead mRead){
    count = (short)mRead.count;
    start = (short)mRead.start;
    end = (short)mRead.end;
  }

  void extendAlignmentLeft(Interval gap){
    byte[] gapSeq = gap.concat;
    short junction = gap.junction;
    end ++;
    if(mappingType == 0){
      // read is aligned to left edge
      if(seq[aln.junctionLeft + 1] == gapSeq[junction]){
        aln.end1 ++;
        aln.leftIdentity ++;
        aln.leftAlnLen ++;
        aln.junctionLeft ++;
      }else{
        aln.junctionLeft = -1;
        aln.leftAlnLen = 0;
        aln.leftIdentity = 0;
      }
    }else if(mappingType == 2){
      // read is aligned on both edges
      if(aln.junctionLeft + 1 == aln.junctionRight){
        //bad alignment -> pick best side
        if(aln.leftIdentity >= aln.rightIdentity){
          if(seq[aln.junctionRight] == gapSeq[junction]){
            aln.junctionLeft ++;
            aln.end1 = aln.junctionRight + 1;
            aln.junctionRight = -1;
            aln.leftAlnLen ++;
            aln.leftIdentity ++;
            aln.rightIdentity = 0;
            aln.rightAlnLen = 0;
            mappingType = 0;
          }else{
            aln.junctionLeft = -1;
            aln.leftAlnLen = 0;
            aln.leftIdentity = 0;
            aln.start1 = aln.junctionRight;
            mappingType = 1;
          }
        }else{
          aln.start1 = aln.junctionRight;
          aln.junctionLeft = -1;
          aln.leftIdentity = 0;
          aln.leftAlnLen = 0;
          mappingType = 1;
        }
      }else{
        // normal case
        if(seq[aln.junctionLeft + 1] == gapSeq[junction]){
          aln.leftIdentity ++;
          aln.junctionLeft ++;
          aln.leftAlnLen ++;
        }else{
          aln.start1 = aln.junctionRight;
          aln.junctionLeft = -1;
          aln.leftIdentity = 0;
          aln.leftAlnLen = 0;
          mappingType = 1;
        }
      }
    }
  }

  void extendAlignmentRight(Interval gap){
    byte[] gapSeq = gap.concat;
    short junction = gap.junction;
    if(mappingType == 1){
      // read is aligned to right edge
      if(seq[aln.junctionRight - 1] ==
          gapSeq[junction + 1]){
        aln.start1 --;
        end ++;
        aln.rightIdentity ++;
        aln.rightAlnLen ++;
        aln.junctionRight --;
      }else{
        aln.junctionRight = -1;
        aln.rightAlnLen = 0;
        aln.rightIdentity = 0;
      }
    }else if(mappingType == 2){
      // read is aligned on both edges
      if(aln.junctionLeft + 1 == aln.junctionRight){
        //bad alignment -> pick best side
        if(aln.rightIdentity >= aln.leftIdentity){
          if(seq[aln.junctionLeft] == gapSeq[junction + 1]){
            aln.start1 = aln.junctionLeft;
            aln.junctionLeft = -1;
            aln.rightAlnLen ++;
            aln.leftIdentity = 0;
            aln.rightIdentity ++;
            aln.junctionRight --;
            end++;
            mappingType = 1;
          }else{
            aln.end1 = aln.junctionRight;
            aln.junctionRight = -1;
            aln.rightIdentity = 0;
            aln.rightAlnLen = 0;
            end ++;
            mappingType = 0;
          }
        }else{
          aln.end1 = aln.junctionRight;
          aln.junctionRight = -1;
          aln.rightIdentity = 0;
          aln.rightAlnLen = 0;
          end ++;
          mappingType = 0;
        }
      }else{
        // normal case
        if(seq[aln.junctionRight - 1] == gapSeq[junction + 1]){
          end ++;
          aln.rightIdentity ++;
          aln.junctionRight --;
          aln.rightAlnLen ++;
        }else{
          aln.end1 = aln.junctionRight;
          aln.junctionRight = -1;
          aln.rightIdentity = 0;
          aln.rightAlnLen = 0;
          end ++;
          mappingType = 0;
        }
      }
    }
  }

  void alignSW(Interval gap, Aligner aligner){
    byte[] g = Arrays.copyOfRange(gap.concat, start, end);
    aln = new EdgeAlignment(aligner.align(seq, g), start, gap);
  }

  void alignSWLeft(Interval gap, Aligner aligner){
    byte[] g = Arrays.copyOfRange(gap.left, start, end);
    aln = new EdgeAlignment(aligner.align(seq, g), start, gap.left.length, true);
  }

  void alignSWRight(Interval gap, Aligner aligner){
    byte[] g = Arrays.copyOfRange(gap.right, start, end);
    aln = new EdgeAlignment(aligner.align(seq, g), start, gap.right.length, false);
  }

}
