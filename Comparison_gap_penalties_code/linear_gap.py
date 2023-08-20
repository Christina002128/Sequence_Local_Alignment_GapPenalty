# linear gap local alignment
#Modified from Simon Tomlinson Bioinformatics Algorithms

def create_matrix(rows, cols):
    my_matrix = [[0 for col in range(cols+1)] for row in range(rows+1)]
    return my_matrix

def calc_score(matrix, x, y):
   sc = seqmatch if sequence1[y- 1] == sequence2[x - 1] else seqmismatch
   base_score = matrix[x - 1][y - 1] + sc
   insert_score = matrix[x - 1][y] + seqgap
   delete_score = matrix[x][y - 1] + seqgap  
   v=max(0, base_score, insert_score, delete_score)
   return v

#makes a single traceback step
def traceback(mymatrix,maxv):
    x=maxv[0]
    y=maxv[-1]
    val=mymatrix[x][y]
    sc = seqmatch if sequence1[y - 2] == sequence2[x - 2] else seqmismatch
    base_score = mymatrix[x - 1][y - 1] + sc
    if base_score==val:
        return [x-1,y-1]

    insert_score = mymatrix[x - 1][y] + seqgap
    if insert_score==val:
        return [x-1,y]
    else:
        return [x,y - 1]


#builds the initial scoring matrix used for traceback
def build_matrix(mymatrix):
    rows=len(mymatrix)
    cols=len(mymatrix[0])
    for i in range(1, rows):
         for j in range(1, cols):
          mymatrix[i][j]  = calc_score(mymatrix, i, j)
    return mymatrix

#gets the max value from the built matrix
def get_max(mymatrix):
     max=mymatrix[0][0]
     mrow=0
     mcol=0
     rows = len(mymatrix)
     cols = len(mymatrix[0])
     for i in range(1, rows):
         for j in range(1, cols):
             if mymatrix[i][j]>max:
                 max=mymatrix[i][j]
                 mrow=i
                 mcol=j
     print("max score: ",max)
     return [mrow,mcol]

#print out the best scoring path from the SW matrix
def print_matrix(mymatrix):
    rows = len(mymatrix)
    cols = len(mymatrix[0])
    s1="  " +sequence1
    s2=" "+sequence2
    print("Dimensions: r= %2d , c= %2d" % (rows, cols))
    for a in s1:
        print(a, end ="")
        print(" \t", end ="")
    print("\n",end="")
    for i in range(0, rows):
        print(s2[i], end ="")
        print(" \t", end ="")
        for j in range(0, cols):
           print("%02d\t" % (mymatrix[i][j]),end="")
        print("\n",end="")

#print out the traceback of the best scoring alignment
def print_traceback(mymatrix):
    #this will print as expected with internal gaps
   print("Building traceback...")
   maxv=get_max(mymatrix)
   print(maxv)
   topstring=""
   midstring=""
   bottomstring=""
   #pad the sequences so indexes into the sequences match the matrix indexes
   asequence1 = "#" + sequence1
   asequence2 = "#" + sequence2
   topstring += asequence1[maxv[-1]]
   bottomstring += asequence2[maxv[0]]
   if asequence1[maxv[-1]] == asequence2[maxv[0]]:
       midstring += "|"
   else:
       midstring += " "
   old_maxv=maxv

   #add the rest of the elements
   search=True
   while(search):
       maxv=traceback(mymatrix,maxv)
       if(mymatrix[maxv[0]][maxv[-1]]==0):
           search=False
           continue

        #deal with single gaps or matches
       if(old_maxv[-1]==maxv[-1]):
           topstring+="-"
           bottomstring +=asequence2[maxv[0]]
       elif(old_maxv[0] == maxv[0]):
            # insertion or deletion
           topstring += asequence1[maxv[-1]]
           bottomstring += "-"
       else:
           # add the next element and go to previous
           topstring += asequence1[maxv[-1]]
           bottomstring += asequence2[maxv[0]]

       if (asequence1[maxv[-1] ] == asequence2[maxv[0] ]) & (old_maxv[-1]!=maxv[-1]) & (old_maxv[0] != maxv[0]):
           midstring += "|"
       else:
        if bottomstring[-1]=='-' or topstring[-1]=='-':
            midstring+=" "
        else :
            midstring += "."
       old_maxv = maxv

   print(topstring[::-1])
   print(midstring[::-1])
   print(bottomstring[::-1])
   return midstring


#build the SW alignment...
def perform_smith_waterman(seq1,seq2):
  global  seqmatch
  global  seqmismatch
  global  seqgap
  global sequence1
  global sequence2
  seqmatch =10
  seqmismatch=-5
  seqgap=-10
#   sequence1 = "GTGTATTTTTTT"
#   sequence2 = "AAAAGTGTTATT"
#input sequences
  sequence1=seq1
  sequence2=seq2
  print("Sequence1: "+sequence1);
  print("Sequence2: "+sequence2);
  mymatrix=create_matrix(len(sequence2), len(sequence1))
  mymatrix=build_matrix(mymatrix)
#   print_matrix(mymatrix)
  midstring=print_traceback(mymatrix)
  total_count=len(midstring)
  mismatch_count=midstring.count('.')
  gap_count=midstring.count(' ')
  match_count=midstring.count('|')
  print(f'alignment length: {total_count}, match:{match_count}, mismatch:{mismatch_count}, gap:{gap_count}')

def read_fasta_filename(filename):
   seq=""
   with open(filename, 'r') as filehandle:
       for line in filehandle:
           if line[0]=='>':
            continue
           seq=seq+line.replace('"','')
       return seq.replace('\n', '')

def argsparse():
    import argparse # command program_name seq1file seq2file int_number
    parser = argparse.ArgumentParser(description='Aligning sequences...') # creating an ArgumentParser object
    parser.add_argument('seq1',action="store",help="First sequence")
    parser.add_argument('seq2',action="store",help="Second sequence")
    parser.add_argument('anum',action="store",help="A number",type =int)
    args = parser.parse_args()
    seq1=read_fasta_filename(args.seq1)
    seq2=read_fasta_filename(args.seq2)
    print(seq1,seq2)

    perform_smith_waterman(seq1,seq2)

argsparse()