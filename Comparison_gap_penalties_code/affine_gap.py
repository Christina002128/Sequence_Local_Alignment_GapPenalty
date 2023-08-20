# affine gap local alignment
# Smith Waterman Alignment with implement of affine gaps
# Sourced from Simon Tomlinson Bioinformatics Algorithms Class code except for the extended codes with specific comments
#contains rows lists each of length cols initially set to 0
#index as my_matrix[1][2] my_matrix[R][C]
def create_matrix(rows, cols):
    my_matrix = [[0 for col in range(cols)] for row in range(rows)]
    return my_matrix

#x is row index, y is column index
# follows[r][c]
def calc_score(M_mat,X_mat,Y_mat,T_mat, x, y):
   #print("seq1:",sequence1[y- 1]," seq2: "+sequence2[x - 1],"x:",x," y:",y)
   sc = seqmatch if sequence1[y] == sequence2[x] else seqmismatch
   M_mat[x][y] = max(0, sc + max(M_mat[x-1][y-1],X_mat[x-1][y-1],Y_mat[x-1][y-1]))
   X_mat[x][y]=max(0, opengap+extendgap+M_mat[x][y-1],extendgap+X_mat[x][y-1])
   Y_mat[x][y]=max(0, opengap+extendgap+M_mat[x-1][y],extendgap+Y_mat[x-1][y])
   # a new matix to track the direction, for better traceback
   slist = [M_mat[x-1][y-1],X_mat[x-1][y-1],Y_mat[x-1][y-1]]
   T_mat[x][y] = slist.index(max(slist))

# extended code: modified traceback() function
#makes a single traceback step
def traceback(M_mat,X_mat,Y_mat,T_mat,maxv):
    x=maxv[0]
    y=maxv[1]
    val=maxv[2]
    if T_mat[x][y]==0:
        return [x-1,y-1,M_mat[x-1][y-1]]
    if T_mat[x][y]==1:
        return [x,y-1,M_mat[x][y-1]]
    if T_mat[x][y]==2:
        return [x-1,y,M_mat[x-1][y]]
 
#extended code: revision of build_matrix() function
#builds the initial scoring matrix used for traceback
def build_matrix(rows, cols):
    # three matrices
    M_mat= create_matrix(rows, cols)
    X_mat= create_matrix(rows, cols)
    Y_mat= create_matrix(rows, cols)
    T_mat= create_matrix(rows, cols)
    for i in range(1, rows):
         for j in range(1, cols):
          calc_score(M_mat,X_mat,Y_mat,T_mat, i, j)
    return M_mat,X_mat,Y_mat,T_mat


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
     return [mrow,mcol,max]

#print out the best scoring path from the SW matrix
def print_matrix(mymatrix):
    rows = len(mymatrix)
    cols = len(mymatrix[0])
    print("Dimensions: r= %2d , c= %2d" % (rows, cols))
    print(' ',end='\t')
    for a in sequence1:
        print(a, end ="")
        print(" \t", end ="")
    print("\n",end="")

    for i in range(0, rows):
        print(sequence2[i], end ="")
        print(" \t", end ="")
        for j in range(0, cols):
           print("%02d\t" % (mymatrix[i][j]),end="")
        print("\n",end="")

# extended code: modification of print_traceback() function
#print out the traceback of the best scoring alignment
def print_traceback(M_mat,X_mat,Y_mat,T_mat):
   #this will print as expected with internal gaps
   print("Building traceback...")
   maxv=get_max(M_mat)
   print(maxv)
   #traverse the matrix to find the traceback elements
   #if more than one path just pick one
   topstring=""
   midstring=""
   bottomstring=""

   # add first element and go to previous
   topstring += sequence1[maxv[1]]
   bottomstring += sequence2[maxv[0]]
   if sequence1[maxv[1]] == sequence2[maxv[0]]:
       midstring += "|"
   else:
       midstring += " "

   # here we need to store the old position so we can track if it is an insertion of deletion
   # code assumes highest score cannot be a gap!
   # error: variable equals to another variable might make these two link together !
   old_maxv=[maxv[0],maxv[1],maxv[2]]

   #add the rest of the elements
   search=True
   while(search):
       maxv=traceback(M_mat,X_mat,Y_mat,T_mat,maxv)
       if(0 in maxv):
        search=False
        continue

        #deal with single gaps or matches
       if(old_maxv[1]==maxv[1]):
           topstring+="-"
           bottomstring +=sequence2[maxv[0]]
       elif(old_maxv[0] == maxv[0]):
            # insertion or deletion
           topstring += sequence1[maxv[1]]
           bottomstring += "-"
       else:
           # add the next element and go to previous
           topstring += sequence1[maxv[1]]
           bottomstring += sequence2[maxv[0]]

       if (sequence1[maxv[1] ] == sequence2[maxv[0]]) & (old_maxv[1]!=maxv[1]) & (old_maxv[0] != maxv[0]):
           midstring += "|"
       else:
        if bottomstring[-1]=='-' or topstring[-1]=='-':
            midstring+=" "
        else :
            midstring += "."
       # print(topstring,bottomstring)
       old_maxv = [maxv[0],maxv[1],maxv[2]]
       
# to delete the last mismatch bases
   if midstring != '|':
    topstring=topstring[:-1]
    midstring=midstring[:-1]
    bottomstring=bottomstring[:-1]
   # reverse print
   print(topstring[::-1])
   print(midstring[::-1])
   print(bottomstring[::-1])
   return midstring


#extended code: revision of perform_smith_waterman() function
#build the SW alignment...
def perform_smith_waterman(seq1,seq2):
#values for weights 
# make global variable
  global  seqmatch
  global  seqmismatch
  global  opengap
  global  extendgap
  global sequence1
  global sequence2

#note these are not the exact weights used in the original SW paper
  seqmatch =10
  seqmismatch=-0.5
  # extended code: add extending and opening gap score
  opengap=-10
  extendgap=-0.5
#   sequence1 = "GTGTATTTTTTT"
#   sequence2 = "AAAAGTGTTATT"
#input sequences
  sequence1=seq1
  sequence2=seq2
  print("Sequence1: "+sequence1)
  print("Sequence2: "+sequence2)
  M_mat,X_mat,Y_mat,T_mat=build_matrix(len(sequence2), len(sequence1))
  midstring=print_traceback(M_mat,X_mat,Y_mat,T_mat)
  total_count=len(midstring)
  mismatch_count=midstring.count('.')
  gap_count=midstring.count(' ')
  match_count=midstring.count('|')
  print(f'alignment length: {total_count}, match:{match_count}, mismatch:{mismatch_count}, gap:{gap_count}')
#   print("M matrix:")
#   print_matrix(M_mat)
#   print("X matrix:")
#   print_matrix(X_mat)
#   print("Y matrix:")
#   print_matrix(Y_mat)

def read_fasta_filename(filename):
   seq=""
   with open(filename, 'r') as filehandle:
       for line in filehandle:
           if line[0]=='>':
            continue
           seq=seq+line.replace('"','')
       return seq

def argsparse():
    import argparse # command program_name seq1file seq2file int_number
    parser = argparse.ArgumentParser(description='Aligning sequences...') # creating an ArgumentParser object
    parser.add_argument('seq1',action="store",help="First sequence")
    parser.add_argument('seq2',action="store",help="Second sequence")
    parser.add_argument('anum',action="store",help="A number",type =int)
    args = parser.parse_args()
    seq1=' '+read_fasta_filename(args.seq1).replace('\n','')
    seq2=' '+read_fasta_filename(args.seq2).replace('\n','')
    ##this calls the SW algorithm when the script loads
    perform_smith_waterman(seq1,seq2)

argsparse()