#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define SL 14



typedef struct reads
{
	char *content;
	struct reads* next;
} Read;


/**********************************************/
char reverseCompleteChar(char c)
{
	char ret;
	switch (c)
	{
		case 'A':
			ret = 'T';
			break;
		case 'T':
			ret = 'A';
			break;
		case 'C':
			ret = 'G';
			break;
		case 'G':
			ret = 'C';
			break;
		default:
			ret = 'N';
			break;
	}
	return ret;
}   
/**********************************************/
void reverseComplete (char *seq, char *rcSeq , int length)
{
	int i;
	for (i=0; i<length; i++)
	{
		rcSeq[i]=reverseCompleteChar (seq[length-1-i]) ;
	}
}
/**********************************************/
int hash (char *s, int b, int l)
{
	int j;
	long x=0;
	for (j=0; j<l; j++)
	{
		//s[j]=toupper(s[j]);
		if (s[(j+b) % l] == 'N' || s[(j+b) % l]=='n')
			return -1;
		else 
		{
			switch (s[(j+b) % l])
			{
				case 'a':
				case 'A': x = (x << 2) | 0; break;
				case 'c':
				case 'C': x = (x << 2) | 1; break;
				case 'g':
				case 'G': x = (x << 2) | 2; break;
				case 't':
				case 'T': x = (x << 2) | 3; break;
			}
		}
	}
	return x;

}
/**********************************************/
int rep(int start, int end, int num)
{
	int shift = (13-end)*2;
	int mask = (int) pow(4, end-start+1)-1;
	return ( (mask<<shift)& num ) >> shift;
}


int main(int argc, char* argv[])
{
	int *seeds = (int *) malloc (256*sizeof(int));
	int *cnt = (int *) malloc ((int)pow(4,14)*sizeof(int));
	int *bb = (int *) malloc ((int)pow(4,14)*sizeof(int));
	int i;
	int x,y,z,w,v,u;
	


	for (i=0; i<256; i++)
	{
		//fprintf(stdout, "%d\n", i>>6);
		x = i>>6 & 3;
		y = i>>4 & 3;
		z = i>>2 & 3;
		w = i & 3;
		if ( x < y && y > z )
		{
			//fprintf (stdout, "1 [%d] %d %d %d %d \n", i, x, y, z, w);
			seeds[i] = 1;
		}

		if ( y==z && x!=y && y!=w)
		{
			//fprintf (stdout, "2 [%d] %d %d %d %d \n", i, x, y, z, w);
			seeds[i]=1;
		}

		if (x!=y && y<z && z<w )
		{
			//fprintf (stdout, "3 [%d] %d %d %d %d\n", i, x, y, z, w);
			seeds[i]=1;
		}

		if (x==y && y==z && z!=w)
		{
			//fprintf (stdout, "4 [%d] %d %d %d %d\n", i, x, y, z, w);
			seeds[i]=1;
		}
		
		/*if (x==y && y==z && z==w)
		{
			seeds[i]=1;
		}*/
	}


	int count=0;
	for (i=0; i<(int)pow(4,14); i++)
	{	
		cnt[i]=0;
		//i = 263172;
		int start[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
		int end[10]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
		int pos = 0;
		int marker = 0;
		int markerCnt = 0;
		while(pos < 12)
		{
			x=i >> (13-pos)*2 & 3;
			y=i >> (12-pos)*2 & 3;
			z=i >> (11-pos)*2 & 3;
			w=i >> (10-pos)*2 & 3;
			v=i >> (9-pos)*2 & 3;
			u=i >> (8-pos)*2 & 3;

			//fprintf(stdout, "%d %d] %d %d %d %d %d %d\n",i, pos, x, y, z, w, v, u);
			markerCnt++;
			//RULE 1
			if (pos < 12 &&  x<y && y >z)
			{
				
				start[markerCnt]=pos; end[markerCnt]=pos+2;

				//fprintf(stdout, "R1 %d\n", y);
				pos+=1;
				marker = marker<<2|y;
			}
			else // RULE 2
			if (pos<10 && x!=y && y==z && z==w && w!=v)
			{
				start[markerCnt]=pos; end[markerCnt]=pos+4;
				//fprintf(stdout, "R2 %d\n", y);
				pos+=1;
				marker = marker<<2|y;
			}
			else /// RULE 3
			if (pos<11 && x<y && y==z && y>w)
			{
				start[markerCnt]=pos; end[markerCnt]=pos+3;
				//fprintf(stdout, "R3 %d\n", y);
				pos+=1;
				marker = marker<<2|y;
			}
			else // RULE 4
			if (pos<9 && x!=y && y==z && z==w && w==v && v!=u)
			{
				start[markerCnt]=pos; end[markerCnt+1]=end[markerCnt]=pos+5;
				//fprintf(stdout, "R4 %d\n", y);
				pos+=1;
				markerCnt++;
				marker = marker<<4|(y*4+y);
			}
			else // RULE 5
			if (pos<10 && x>=y && y>z && z<w && w<=v)
			{
				start[markerCnt]=pos; end[markerCnt]=pos+4;
				//fprintf(stdout, "R5 %d\n", y);
				pos+=1;
				marker = marker<<2|z;
			}
			else // RULE 0
			if (pos<10&&  (x==y && y==z && z==w && w==v)) // 
			    //y==z && z==w && w==v && v==u) )
			{
			
				start[markerCnt]=pos; end[markerCnt]=pos+4;
				marker = x*(int)(pow(4,0)+pow(4,1)+pow(4,2)+pow(4,3)+pow(4,4));
				if (marker == 0 ) marker = 1;
				fprintf(stdout, "%d %d %d %d %d %d\n", i, pos, pos+4, marker, 4, rep(pos, pos+4, i) );
				//fprintf(stdout, "R6 %d\n", y);
				//bb[i]=marker+1;
				count++;
				markerCnt=4;
				break;
			}
			else
			{
				pos++;
				markerCnt--;
			}

			//fprintf(stdout, "===> %d %d %d \n", markerCnt, start[markerCnt], end[markerCnt]);
			if (markerCnt >= 4)
			{
				if (markerCnt == 4)
				{
					if (seeds[marker]==1)
					{
						count++;
						fprintf(stdout, "%d %d %d %d %d %d\n", i, start[1], end[4], marker,(end[4]-start[1]+1) ,  rep(start[1], end[4], i) );
						bb[i]=rep(start[1], end[4], i)+1;
						break;
					}
					else 
					{
						if (start[2] == -1)
						{
							for (z=3; z<10; z++)
							{
								start[z-2]=start[z];
								end[z-2]=end[z];
							}
							markerCnt = 2;
							marker = marker & 15;
						}
						else
						{
							for (z=2; z<10; z++)
							{
								start[z-1]=start[z];
								end[z-1]=end[z];
							}

							marker = marker & 63;
							markerCnt=3;
						}
					}
				}
				else if (markerCnt == 5)
				{
					if (seeds[marker >> 2 & 255] == 1)
					{
						count++;
						fprintf(stdout, "%d %d %d %d %d %d\n", i, start[2], end[5], marker, (end[5]-start[2]+1),rep(start[2], end[5], i) );
						bb[i]=rep(start[2], end[5], i)+1;
						break;
					}
					else if (seeds[marker & 255]==1 && start[2]!=-1)
					{
						count++;
						fprintf(stdout, "%d %d %d %d %d %d\n", i, start[1], end[4], marker, (end[4]-start[1]+1), rep(start[1], end[4], i) );
						bb[i]=rep(start[1], end[4], i)+1;
						break;
					}
					else
					{
						if (start[2]==-1)
						{
							for (z=3; z<10; z++)
							{
								start[z-2]=start[z];
								end[z-2]=end[z];
							}
							markerCnt = 3;
							marker = marker & 31;
						}
						else
						{
						
							for (z=4; z<10; z++)
							{
								start[z-3]=start[z];
								end[z-3]=end[z];
							}
							markerCnt = 2;
							marker = marker & 15;
						}
					}
				}
			}
		}
		//break;
	}

//	return 1;
//	fprintf (stdout, "=====> %d <=====", count);


	char *s= (char *) malloc(100);
	char *rs = (char *)malloc(100);
	int h;
	int b;
	int j;
	int f=0;
	int ll;
	int mm=0;

	int lpos, lmax;
	while (fgets(s, 1000, stdin)!=NULL)
	{
		lpos = -1;
		lmax = -1;
		f=0;
		/*ll = strlen(s);
		if (++mm % 2)
		{
			reverseComplete(s, rs, ll);
			s = rs;
			rs[ll]='\n';
		}*/
		for (i=0; i<strlen(s)-14;i++)
		{

			h = hash (s+i, 0, 14);
			if (h>=0 && bb[h])
			{
				/*if (cnt[h]>lmax)
				{
					lpos = h;
					lmax = cnt[h];
				}*/
				f=1;
				break;
			}
		}
		if (f)
		{
			//cnt[lpos]++;
			fprintf(stdout, "%d %s",  bb[h], s);
		}
		else
		{
			//fprintf(stdout, "-1 %s", s);
		}
	}

	/*for (i=0; i<(int)pow(4,14); i++)
	{
		if (cnt[i]>0)
			fprintf(stdout, "%d\n", cnt[i]);
	}*/



}
