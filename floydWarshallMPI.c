#include <mpi.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>


#define SIZE_INIT 10
#define FILE_ERROR 2
#define MEMORY_ERROR 3
#define ROOT 0
#define NEXT(r, p) (r+1)%p
#define PREVIOUS(r, p) (r+p-1)%p
#define FILE_IS_NULL(FILE) if ( FILE == NULL ) { printf("Erreur fichier introuvable\n");MPI_Abort(MPI_COMM_WORLD,FILE_ERROR);MPI_Finalize(); exit(FILE_ERROR); }
#define MIN(x,y) (((x) < (y)) ? (x) : (y))


//suposition multiplication de matrice (n,m)*(m,k) = (n,k) => nombre de ligne a circuler = nombre de colonnes a circuler ..ouf

typedef struct{
      int column;
      int line;
      int *tab;
    }Matrix;
    
typedef struct{
  Matrix* mat;
  int index;
}IndexedMatrix;
  
void printMatrix(Matrix* mat);

typedef struct Node {
    IndexedMatrix* item;
    struct Node* next;
} Node;

typedef struct File {
    Node* head;
    Node* tail;
    int size;
} File;

void addElement(File* file , IndexedMatrix* elem){
  
  Node* n =  malloc (sizeof(Node));
    n->item = elem;
    n->next = NULL;

    if (file->head == NULL) { // no head
        file->head = n;
	file->tail=n;
    } else{
        file->tail->next = n;
    }
    file->tail = n;
    file->size++;
}
IndexedMatrix* getElemAt(File* file, int index){
  if(index> file->size){
    printf("error too big index\n");
    MPI_Abort(MPI_COMM_WORLD,FILE_ERROR);
  }
  int i=0;
  Node* head = file->head;
  while(i<index){
    head=head->next;
    i++;
  }
   return head->item;
}

IndexedMatrix* popElement(File* queue) {
    // get the first item
    Node* head = queue->head;
    IndexedMatrix* item = head->item;
    queue->head = head->next;
    queue->size--;
    // free the memory of original head
    free(head);
    return item;
}

int sizeOf(File* file){
  return file->size;
}
int world_rank,world_size,previous,next,value,lineOnRequest;
MPI_Group custom,origin;
MPI_Comm new_comm;



Matrix* allocateMatrix(int line,int column ){
  Matrix* mat = malloc(sizeof(Matrix));
  if(mat==NULL){
    exit(MEMORY_ERROR);
  }
  mat->tab=malloc((line*column)*(sizeof(int)));
  if(mat->tab==NULL){
    exit(MEMORY_ERROR);
  }
  mat->column=column;
  mat->line=line;
  int elem;
  for(elem=0;elem<(mat->line)*(mat->column);elem++){
	mat->tab[elem]=INT_MAX;
  }
  return mat;
}
void printMatrix(Matrix* mat){
  int elem;
  for(elem=0;elem<(mat->line)*(mat->column);elem++){
    if(elem%(mat->column)==0 && elem!=0){
      printf("\n");
    }
    printf(" %d ",mat->tab[elem]);
  }
   printf("\n");
}

int getSizeOfLine(char* nameFile){
  int size=0,value;
  char* charElement = NULL;
  size_t len = 0;
  FILE* myfile=fopen(nameFile, "r");
  FILE_IS_NULL(myfile)
  getline(&charElement, &len, myfile);
  int i=0;
  while (sscanf(charElement, "%d%n ", &value,&i)==1) {
    charElement+=i;
    size++;
  }
      fclose(myfile);
      return size;
}

void allocateNewLine(Matrix* mat,int cmpt){
  mat->tab =  realloc(mat->tab,(cmpt + mat->column)*sizeof(int));
  if(mat->tab==NULL){
    exit(MEMORY_ERROR);
  }
  mat->line++;
  while(cmpt<(mat->line)*(mat->column)-1){
    mat->tab[cmpt++]= -1;
  }
}

Matrix* readMatrix(char* nameFile){
    int value,cmpt=0;
    Matrix* mat = malloc(sizeof(Matrix));
    mat->column=getSizeOfLine(nameFile);
    mat->line=1;
    FILE* myfile; 
    myfile = fopen(nameFile, "r");
    FILE_IS_NULL(myfile)
    mat->tab= malloc(sizeof(int)*mat->column);
    int diag=0;
    while (fscanf(myfile, "%d", &value)>=1) {
      if((cmpt%(mat->column)) == 0&& (cmpt!=0) ){
	allocateNewLine(mat,cmpt);
	diag++;
      }
      if((cmpt!=(mat->column*diag+diag) && value==0 )){
	mat->tab[cmpt++]=INT_MAX;
      }else{
	mat->tab[cmpt++]=value;
      } 
    }
    fclose(myfile);
    return mat;
}

Matrix* transposeMatrix(Matrix* mat){
  Matrix* Tmat = allocateMatrix(mat->column,mat->line);
  int columnElem=0,lineElem=0;
  for(lineElem=0;lineElem<mat->line;lineElem++){
      for(columnElem=0;columnElem<mat->column;columnElem++){
	Tmat->tab[lineElem*mat->column + columnElem] = mat->tab[columnElem*mat->column + lineElem];
      }
  }
  return Tmat;
}
    
IndexedMatrix* allocateIndexedMatrix(int line ,int  index,Matrix* mat){
  IndexedMatrix* indexMat = malloc(sizeof(IndexedMatrix));
  if(mat==NULL){
    exit(MEMORY_ERROR);
  }
  indexMat->index=index;
  indexMat->mat=allocateMatrix(1,mat->column);
  int i=0;
  while(i<mat->column){indexMat->mat->tab[i]=mat->tab[i+line*mat->column];i++;}
  return indexMat; 
}

int calculation(Matrix* a, Matrix* b){
  int res=INT_MAX;
  int i;
      for(i=0;i<(b->column);i++){
	int vv =(a->tab[i])+(b->tab[i]);
	if((a->tab[i]==INT_MAX) || (b->tab[i]==INT_MAX) ||(b->tab[i]+(a->tab[i])<0)){
	  vv=INT_MAX;
	}
	res = MIN(res,vv);
  }
  return res;
  
}

void createGroupCommunitation(int line){
  MPI_Comm_group(MPI_COMM_WORLD, &origin);
  int arrayCustom[world_size];
  int i=0;
  int lineGap =ceil((double)line/world_size);
  while(line>0){
    if(line<=lineGap){
      lineOnRequest=line;
    }
    line -=lineGap;
    arrayCustom[i]=i;
    i++;
    
  }
  value=i ;
   previous= PREVIOUS(world_rank,i);
   
  next=NEXT(world_rank,i);
  MPI_Group_incl(origin, value, arrayCustom, &custom);
  MPI_Comm_create(MPI_COMM_WORLD, custom, &new_comm); 
}


File* initFile(){
  File * k=malloc(sizeof(File));
  k->head=NULL;
  k->tail=NULL;
  k->size=0;
  return k;
}

File* putInFile(int* flagUsingLine,Matrix* mat){
  int line=0;
  File* file= initFile();
//   file->elem = malloc(mat->line * sizeof(IndexedMatrix*));
  if(file==NULL){MPI_Abort(MPI_COMM_WORLD,MEMORY_ERROR);MPI_Finalize(); exit(MEMORY_ERROR);}
  for(line=0; line < mat->line;line++){
    addElement(file, allocateIndexedMatrix(line , flagUsingLine[line], mat));
  }
  return file;
}

void Send(IndexedMatrix* matIndex,int next){
  MPI_Send(&(matIndex->mat->line), 1, MPI_INT, next, 0, new_comm);
  MPI_Send(&(matIndex->mat->column), 1, MPI_INT, next, 1, new_comm);
  MPI_Send(matIndex->mat->tab, (matIndex->mat->column)*(matIndex->mat->line), MPI_INT, next, 2, new_comm);
  MPI_Send(&(matIndex->index), 1, MPI_INT, next, 3, new_comm);
}

IndexedMatrix* Receiv(int previous){
  
  int line , column , index;
  IndexedMatrix* iMat;
  MPI_Recv(&line , 1,MPI_INT,previous,0,new_comm,MPI_STATUS_IGNORE);
  MPI_Recv(&column , 1,MPI_INT,previous,1,new_comm,MPI_STATUS_IGNORE);
  Matrix* k= allocateMatrix(line,column);
  MPI_Recv(k->tab , line*column,MPI_INT,previous,2,new_comm,MPI_STATUS_IGNORE);
  MPI_Recv(&index , 1,MPI_INT,previous,3,new_comm,MPI_STATUS_IGNORE);
  iMat=allocateIndexedMatrix(0,index,k);
  return iMat;
}

Matrix* cpy(Matrix * mat){
  Matrix * cpy = allocateMatrix(mat->line, mat->column);
  int k=0;while(k<mat->line*mat->column){cpy->tab[k]=mat->tab[k];k++;}
  return cpy;
} 

void formatMatrix(Matrix* mat){
  int elem;
  #pragma omp parallel for
  for(elem=0;elem < mat->line*mat->column;elem++) {
    if(((elem/mat->column)!=elem%(mat->line))&& mat->tab[elem]==INT_MAX){
      mat->tab[elem]=-1;
    }
  }
}


int main(int argc, char **argv) {
  
  int totColumn=0,totLine=0,wild=0;
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  File* fileMat;
  File* fileTMat;
  if(world_size==1){
    printf("C'est pas un code séquentiel M^ù! !!!!");
    MPI_Abort(MPI_COMM_WORLD,2);
    MPI_Finalize();
  }
  Matrix* mat;
  Matrix* Tmat;
 
  
  if(world_rank==ROOT){
    mat =  readMatrix(argv[1]);
    Tmat = transposeMatrix(mat);
    totColumn=Tmat->column;
    totLine= mat->line;
    wild = mat->column;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&totLine, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&totColumn, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&wild, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    
  int lineFlags[totLine];// Represent at each time what lines are present (mat)
  int columnFlags[totColumn];// Represent at each time what lines are present (Tmat)
  if(world_rank==ROOT){
    int flag=0;
    for(flag=0;flag<totLine && flag< totColumn; flag++){
      if(totColumn>flag)columnFlags[flag]=flag;
      if(totLine>flag)lineFlags[flag]=flag;
    }
  }
  //<=> total colonne car (n,m)*(m,k) ligne matrice = colonnes Tmatrice
  //Reduit le groupe de communication si nb ligne < nb proc
  createGroupCommunitation(totLine);
//On ne parle qu'au groupe réduit
  if(world_rank < value){
   int flagUsingLine[(int)ceil((double)totLine/value)];
   int flagUsingColumn[(int)ceil((double)totLine/value)];
    
    
    if(world_rank!= ROOT){
      int nbLigne=ceil((double)totLine/value);
      int nbColumnTmat=ceil((double)wild/value);
      mat=allocateMatrix(nbLigne,totColumn);
      Tmat=allocateMatrix(nbColumnTmat,totColumn);
    }
    else{
      int elem=0;
      while(elem < ceil((double)totLine/value)){
	  flagUsingLine[elem]=elem;
	  flagUsingColumn[elem]=elem;
	elem++;
      }
    }
    int currentPow=0;
      File* power= initFile();
    IndexedMatrix* indexedmatInit = malloc(sizeof(IndexedMatrix));
    indexedmatInit->index=1;
    indexedmatInit->mat = mat;
    if(world_rank==ROOT){
      addElement(power,indexedmatInit);
    }
    while(currentPow!=totLine){
      printf("currentPow %d \n",currentPow);
      if(world_rank==ROOT){
	if(currentPow!=0){
	  int indice,sizepower =sizeOf(power),matPow=getElemAt(power,sizepower-1)->index;
	  int matToPow=getElemAt(power,sizepower-1)->index;
	  if(matPow+matToPow>totLine){
	     int k;
	     for(k=0;k<sizepower;k++){
	      if(currentPow+ getElemAt(power,k)->index<=totLine){
		currentPow=matPow+getElemAt(power,k)->index;
		indice=k;
	      }
	    }
	  }else{ 
	    indice=sizepower-1;
	    currentPow=2*matPow;
	  }
	  mat=getElemAt(power,sizeOf(power)-1)->mat;
	  Tmat= transposeMatrix(getElemAt(power,indice)->mat);
	}
	else{
	  currentPow=2;
	}
      }
      //Distribution des donnees
      MPI_Scatter(mat->tab, ceil((double)totLine/value)*totColumn, MPI_INT, mat->tab,ceil((double)totLine/value)*totColumn, MPI_INT, ROOT, new_comm);
      MPI_Scatter(Tmat->tab, ceil((double)wild/value)*totColumn, MPI_INT,Tmat->tab,ceil((double)wild/value)*totColumn, MPI_INT, ROOT, new_comm); 
      
      //Distribution des flags
	int taille=ceil((double)totLine/value);
	if(world_rank==ROOT){
	  int rank=0,indice;
	  for(rank=1;rank<value;rank++){
	    indice = taille+taille*(rank-1);//debut partie tableau pour ligne/colonne a distribuer
	    int tmpTaille = taille;
	    if(rank==value-1){
	      tmpTaille= lineOnRequest;
	  }
	  MPI_Send(lineFlags+indice, tmpTaille, MPI_INT, rank, 0, new_comm);
	  MPI_Send(columnFlags+indice, tmpTaille, MPI_INT, rank, 1, new_comm);
	}
	
	//Maj de ses données
	mat->line=taille;
	Tmat->line=taille;
	int k=0;
	while(k<taille){flagUsingLine[k]=k;k++;}
	
	
      }else{
	int tmpTaille = taille;
	  if(world_rank==value-1){
	    tmpTaille= lineOnRequest;
	  }
	MPI_Recv(&flagUsingLine , tmpTaille,MPI_INT,ROOT,0,new_comm,MPI_STATUS_IGNORE);
	MPI_Recv(&flagUsingColumn , tmpTaille,MPI_INT,ROOT,1,new_comm,MPI_STATUS_IGNORE);
	
	if(world_rank==value-1){
	  mat->line = lineOnRequest;
	  Tmat->line=lineOnRequest;
	}
	
      }
	//initiation file
	fileMat=putInFile(flagUsingLine,mat);
	fileTMat=putInFile(flagUsingColumn,Tmat);
	
	//boucle cal+ circulation
	int res [totLine*totColumn];
	int LKIN;
	for(LKIN=0;LKIN<totLine*totColumn;LKIN++){res[LKIN]=0;}
	 
	int circul;
	for(circul=0;circul<totLine;circul++){

	  int stackElem;
	  int size=sizeOf(fileMat);
	  #pragma omp parallel for
	  for(stackElem=0;stackElem<size;stackElem++){
	    IndexedMatrix* elemFileMat = getElemAt(fileMat,stackElem);
	    IndexedMatrix* elemFileTMat = getElemAt(fileTMat,stackElem);
	    Matrix * transposed = transposeMatrix(elemFileTMat->mat);
	    res[elemFileMat->index*totColumn+elemFileTMat->index]= calculation(elemFileMat->mat,elemFileTMat->mat);;
	    free(transposed);
	    }
	    IndexedMatrix* mat=popElement(fileTMat);
	     IndexedMatrix* mt;
	     if(world_rank%2 ==0){
	       Send(mat,next);
	       mt=Receiv(previous);
	    }
	    else{
	      mt=Receiv(previous);
	      Send(mat,next);
	    }
	     addElement(fileTMat,mt);
	    MPI_Barrier(new_comm);
	}
	//Gather Reduce
	int* tmp = malloc(sizeof(int)*totLine*totColumn);
	if(tmp==NULL){MPI_Abort(MPI_COMM_WORLD,MEMORY_ERROR);}
	 MPI_Reduce( res, tmp, totLine*totColumn, MPI_INT, MPI_SUM, ROOT, new_comm);
	 if(world_rank==ROOT){
	 free(mat->tab);
	 mat->tab=tmp;
	 mat->line=totLine;
	 IndexedMatrix* i= malloc(sizeof(IndexedMatrix));
	 i->index=currentPow;
	 i->mat=cpy(mat);
	 addElement(power,i);
	}
	MPI_Bcast(&currentPow, 1, MPI_INT, ROOT, new_comm);
      }
      printf("currentPow %d \n",currentPow);
      if(world_rank==ROOT){
	formatMatrix(mat);
	printMatrix(mat);
      }
    }
  MPI_Finalize();
  return 0;
}