/**
 * Created by Геннадий on 05.08.2014.
 */
class EdgeAlignment extends Constants{

  int start1;
  int end1;

  int junctionLeft;
  int junctionRight;
  int leftIdentity;
  int rightIdentity;
  int leftAlnLen;
  int rightAlnLen;

  EdgeAlignment(Alignment aln, int readStart, Interval gap){
    int junction = gap.junction;
    byte[] concatSeq = gap.concat;
    start1 = aln.start1;
    end1 = aln.end1;
    int alnPos = 0;
    int readPos = start1;
    int genomePos = readStart + aln.start2;
    leftIdentity = 0;
    rightIdentity = 0;
    int rightAlnPos = 0;
    if(readStart + aln.end2 <= junction || genomePos > junction + 1){
      // alignment is useless
      junctionLeft = -1;
      junctionRight = -1;
      rightAlnLen = 0;
      leftAlnLen = 0;
      return;
    }
    if(genomePos <= junction){
      for(; alnPos < aln.length && genomePos <= junction; alnPos++){
        byte c = aln.sequence1[alnPos];
        if(c == GAP){
          genomePos++;
        }else if(c < 'a'){
          if(c == concatSeq[genomePos]){
            leftIdentity++;
          }
          readPos++;
          genomePos++;
        }else{
          readPos++;
        }
      }
      if(aln.sequence1[alnPos - 1] == GAP){
        junctionLeft = -1;
        leftAlnLen = 0;
        start1 = readPos;
      }else{
        junctionLeft = readPos - 1;
        leftAlnLen = alnPos;
      }
    }else{
      junctionLeft = -1;
    }
    if(alnPos == aln.length){
      junctionRight = -1;
      rightAlnLen = 0;
    }else{
      for(; alnPos < aln.length && genomePos <= junction + 1; alnPos++){
        byte c = aln.sequence1[alnPos];
        if(c == GAP){
          genomePos++;
        }else if(c < 'a'){
          if(c == concatSeq[genomePos]){
            rightIdentity ++;
          }
          readPos ++;
          genomePos ++;
        }else{
          readPos ++;
        }
      }
      if(aln.sequence1[alnPos - 1] == GAP){
        rightAlnLen = 0;
        junctionRight = -1;
        end1 = readPos;
      }else{
        junctionRight = readPos - 1;
        rightAlnLen = aln.length - alnPos + 1;
        rightAlnPos = alnPos - 1;
      }
      for(; alnPos < aln.length; alnPos++){
        byte c = aln.sequence1[alnPos];
        if(c == GAP){
          genomePos++;
        }else if(c < 'a'){
          if(c == concatSeq[genomePos]){
            rightIdentity ++;
          }
          readPos ++;
          genomePos ++;
        }else{
          readPos ++;
        }
      }
    }
    if(junctionLeft + 1 == junctionRight){
      //look for insertion right to junction
      boolean isInsert = false;
      int insertionStartAlnPos = -1;
      for(alnPos = rightAlnPos; alnPos < aln.length; alnPos++){
        byte c = aln.sequence1[alnPos];
        if(c != GAP){
          if(c < 'a'){
            if(isInsert){
              //insertion end
              int insertionEnd = alnPos;
              int insertionLen = insertionEnd - insertionStartAlnPos;
              if((alnPos - rightAlnPos)%insertionLen == 0){
                // tandem repeat check
                int repeatNum = (alnPos - rightAlnPos)/insertionLen;
                boolean equal = true;
                for(int k = 0; k < repeatNum; k++){
                  for(int i = rightAlnPos + k*insertionLen, j = insertionStartAlnPos; j < insertionEnd; i++, j++){
                    if(aln.sequence1[i] != upper[aln.sequence1[j]]){
                      equal = false;
                      break;
                    }
                  }
                }
                if(equal){
                  // move insert left to fill the gap
                  junctionRight += insertionLen;
                  rightAlnLen -= insertionLen;
                }
              }
              // non tandem situation
              if(insertionLen >= insertionStartAlnPos - rightAlnPos){
                boolean equal = true;
                for(int i = insertionStartAlnPos - 1, j = insertionEnd - 1; i >= rightAlnPos; i--, j--){
                  if(aln.sequence1[i] != upper[aln.sequence1[j]]){
                    equal = false;
                    break;
                  }
                }
                if(equal){
                  // move insert left to fill the gap
                  junctionRight += insertionLen;
                  rightAlnLen -= insertionLen;
                }
              }
              break;
            }
          }else{
            if(!isInsert){
              isInsert = true;
              insertionStartAlnPos = alnPos;
            }
          }
        }
      }
      //look for insertion left to junction
      isInsert = false;
      int insertionEndAlnPos = -1;
      for(alnPos = rightAlnPos - 1; alnPos >= 0; alnPos--){
        byte c = aln.sequence1[alnPos];
        if(c != GAP){
          if(c < 'a'){
            if(isInsert){
              //insertion end
              int insertionStart = alnPos;
              int insertionLen = insertionEndAlnPos - insertionStart;
              if((rightAlnPos - alnPos)%insertionLen == 0){
                // tandem repeat check
                int repeatNum = (rightAlnPos - alnPos)/insertionLen;
                boolean equal = true;
                for(int k = 0; k < repeatNum; k++){
                  for(int i = rightAlnPos - k*insertionLen, j = insertionEndAlnPos; j >= insertionStart; i--, j--){
                    if(aln.sequence1[i] != upper[aln.sequence1[j]]){
                      equal = false;
                      break;
                    }
                  }
                }
                if(equal){
                  // move insert left to fill the gap
                  junctionLeft -= insertionLen;
                  leftAlnLen -= insertionLen;
                }
              }
              // non tandem situation
              if(insertionLen >= insertionStartAlnPos - rightAlnPos){
                boolean equal = true;
                for(int i = insertionStart, j = insertionEndAlnPos; i >= rightAlnPos; i++, j++){
                  if(upper[aln.sequence1[i]] != aln.sequence1[j]){
                    equal = false;
                    break;
                  }
                }
                if(equal){
                  // move insert left to fill the gap
                  junctionLeft -= insertionLen;
                  leftAlnLen -= insertionLen;
                }
              }
              break;
            }
          }else{
            if(!isInsert){
              isInsert = true;
              insertionEndAlnPos = alnPos;
            }
          }
        }
      }
    }
  }

  EdgeAlignment(Alignment aln, int readStart, int genomeLen, boolean left){
    start1 = aln.start1;
    end1 = aln.end1;
    leftIdentity = 0;
    rightIdentity = 0;
    if(left){
      if(readStart + aln.end2 != genomeLen){
        // alignment is useless
        junctionLeft = -1;
        junctionRight = -1;
        rightAlnLen = 0;
        leftAlnLen = 0;
      }else{
        leftIdentity = aln.identity;
        junctionLeft = end1 - 1;
        junctionRight = -1;
        leftAlnLen = aln.length;
      }
    }else{
      if(readStart + aln.start2 != 0){
        // alignment is useless
        junctionLeft = -1;
        junctionRight = -1;
        rightAlnLen = 0;
        leftAlnLen = 0;
      }else{
        rightIdentity = aln.identity;
        junctionRight = start1;
        junctionLeft = -1;
        rightAlnLen = aln.length;
      }
    }
  }

}
