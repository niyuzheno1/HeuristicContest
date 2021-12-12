#pragma GCC target("avx")
#pragma GCC optimize("O3")
#include <bits/stdc++.h>
using namespace std;

unsigned long xor128() {
  static unsigned long x=123456789, y=362436069, z=521288629, w=88675123;
  unsigned long t=(x^(x<<11));
  x=y; y=z; z=w;
  return (w=(w^(w>>19))^(t^(t>>8)));
}

struct rect{
	int x1,y1,x2,y2;
	int area(){
		return (x2-x1)*(y2-y1);
	}
};

int N,X[200],Y[200],R[200],ord[200],resA[200],resB[200],resC[200],resD[200],vx[10001][200],vy[10001][200],cx[10001],cy[10001],T[200],vra[100][100][200],cra[100][100];

pair<int,int>pv[200];
rect S[200];
double score,maxscore;
constexpr double C0=50000000,C1=100000000000;
constexpr int M=7;

inline bool include(rect r,int p,int q){
	return r.x1<p&&p<r.x2&&r.y1<q&&q<r.y2;
}

inline bool intersect(rect r1,rect r2){
	return (!(r1.x2<=r2.x1||r2.x2<=r1.x1))&&(!(r1.y2<=r2.y1||r2.y2<=r1.y1));
}

inline rect newrect(rect r,int x,int y,int p,int q){
	rect nr;
	if(p<=x){
		nr.x1=p;
		nr.x2=r.x2;
	}
	else{
		nr.x1=r.x1;
		nr.x2=p;
	}
	if(q<=y){
		nr.y1=q;
		nr.y2=r.y2;
	}
	else{
		nr.y1=r.y1;
		nr.y2=q;
	}
	return nr;
}

inline double calc_score(rect &r,int s){
	return 1.0-(1.0-(double)min(r.area(),s)/max(r.area(),s))*(1.0-(double)min(r.area(),s)/max(r.area(),s));
}

inline void addx(int x,int num){
	vx[x][cx[x]]=num;
	cx[x]++;
}

inline void erasex(int x,int num){
	cx[x]--;
	for(int i=0;;i++){
		if(vx[x][i]==num){
			vx[x][i]=vx[x][cx[x]];
			break;
		}
	}
}

inline void addy(int y,int num){
	vy[y][cy[y]]=num;
	cy[y]++;
}

inline void erasey(int y,int num){
	cy[y]--;
	for(int i=0;;i++){
		if(vy[y][i]==num){
			vy[y][i]=vy[y][cy[y]];
			break;
		}
	}
}

inline void rect_update(rect &old_rect,rect &new_rect,int id){
	if(old_rect.x1!=new_rect.x1){
		erasex(old_rect.x1,id);
		addx(new_rect.x1,id);
	}
	if(old_rect.y1!=new_rect.y1){
		erasey(old_rect.y1,id);
		addy(new_rect.y1,id);
	}
	if(old_rect.x2!=new_rect.x2){
		erasex(old_rect.x2,id);
		addx(new_rect.x2,id);
	}
	if(old_rect.y2!=new_rect.y2){
		erasey(old_rect.y2,id);
		addy(new_rect.y2,id);
	}
	old_rect=new_rect;
}

int main(){
	scanf("%d",&N);
	for(int i=0;i<N;i++){
		scanf("%d%d%d",&X[i],&Y[i],&R[i]);
	}
	score=maxscore=0;
	for(int i=0;i<N;i++){
		S[i].x1=X[i];
		S[i].y1=Y[i];
		S[i].x2=X[i]+1;
		S[i].y2=Y[i]+1;
		addx(X[i],i);
		addy(Y[i],i);
		addx(X[i]+1,i);
		addy(Y[i]+1,i);
		pv[i]={R[i],i};
	}
	sort(pv,pv+N);
	for(int i=0;i<N;i++){
		ord[i]=pv[i].second;
	}
	for(int i=0;i<100;i++){
		for(int j=0;j<100;j++){
			for(int k=0;k<N;k++){
				int t=ord[k];
				rect r={min(X[t],(i+1)*100),min(Y[t],(j+1)*100),max(X[t]+1,i*100),max(Y[t]+1,j*100)};
				bool f=true;
				for(int l=0;l<N;l++){
					if(t==l)continue;
					if(intersect(r,S[l])){
						f=false;
						break;
					}
				}
				if(f){
					vra[i][j][cra[i][j]]=t;
					cra[i][j]++;
				}
			}
		}
	}
	for(int _=35000000/pow(N,0.3);_>=0;_--){
		int p,q;
		p=xor128()%10001;
		q=xor128()%10001;
		int pa,qa;
		pa=p/100-(p==10000);
		qa=q/100-(q==10000);
		bool ff=true;
		
		for(int i=0;i<cra[pa][qa];i++){
			int t=vra[pa][qa][i];
			if(include(S[t],p,q)){
				rect nr=newrect(S[t],X[t],Y[t],p,q);
				if(calc_score(nr,R[t])-calc_score(S[t],R[t])>=-(double)_/C0){
					score+=calc_score(nr,R[t])-calc_score(S[t],R[t]);
					rect_update(S[t],nr,t);
				}
				ff=false;
				break;
			}
		}
		if(ff){
			for(int i=0;i<cra[pa][qa];i++){
				int t=vra[pa][qa][i];
				rect nr=newrect(S[t],X[t],Y[t],p,q);
				if(calc_score(nr,R[t])-calc_score(S[t],R[t])>=-(double)_/C0){
					bool f=true;
					for(int j=0;j<N;j++){
						if(j==t)continue;
						if(intersect(nr,S[j])){
							f=false;
							break;
						}
					}
					if(f){
						score+=calc_score(nr,R[t])-calc_score(S[t],R[t]);
						rect_update(S[t],nr,t);
					}
				}
			}
		}
		if(score>maxscore){
			maxscore=score;
			for(int i=0;i<N;i++){
				resA[i]=S[i].x1;
				resB[i]=S[i].y1;
				resC[i]=S[i].x2;
				resD[i]=S[i].y2;
			}
		}
		if(_%4==0){
			p=xor128()%N;
			//-x
			if(S[p].x1!=0){
				rect nr=S[p];
				nr.x1--;
				if(_%M==0&&S[p].y1!=Y[p]){
					nr.y1++;
				}
				if(_%M==1&&S[p].y2!=Y[p]+1){
					nr.y2--;
				}
				int c=0;
				for(int i=0;i<cx[S[p].x1];i++){
					int k=vx[S[p].x1][i];
					if(k==p)continue;
					if(intersect(nr,S[k])){
						T[c]=k;
						c++;
					}
				}
				double d;
				d=calc_score(nr,R[p])-calc_score(S[p],R[p]);
				for(int i=0;i<c;i++){
					rect tmprec=S[T[i]];
					tmprec.x2--;
					if(tmprec.x2<=X[T[i]]){
                      	d-=334334334;
                    }
					else{
						d+=calc_score(tmprec,R[T[i]])-calc_score(S[T[i]],R[T[i]]);
					}
				}
				if(c==0&&d<0&&nr.x2!=X[p]+1){
					nr.x2--;
					while(true){
						if(nr.x1==0||nr.x2==X[p]+1)break;
						bool f=true;
						nr.x1--;
						nr.x2--;
						for(int i=0;i<cx[nr.x1+1];i++){
							if(vx[nr.x1+1][i]==p)continue;
							if(intersect(nr,S[vx[nr.x1+1][i]])){
								f=false;
                              	break;
							}
						}
						if(!f){
							nr.x1++;
							nr.x2++;
							break;
						}
					}
					rect_update(S[p],nr,p);
				}
				else if(d>=-(double)_/C1){
					score+=d;
					rect_update(S[p],nr,p);
					for(int i=0;i<c;i++){
						rect tmprec=S[T[i]];
						tmprec.x2--;
						rect_update(S[T[i]],tmprec,T[i]);
					}
				}
			}
			if(score>maxscore){
				maxscore=score;
				for(int i=0;i<N;i++){
					resA[i]=S[i].x1;
					resB[i]=S[i].y1;
					resC[i]=S[i].x2;
					resD[i]=S[i].y2;
				}
			}
		}		
		if(_%4==1){
			p=xor128()%N;
			//-y
			if(S[p].y1!=0){
				rect nr=S[p];
				nr.y1--;
				if(_%M==0&&S[p].x1!=X[p]){
					nr.x1++;
				}
				if(_%M==1&&S[p].x2!=X[p]+1){
					nr.x2--;
				}
				int c=0;
				for(int i=0;i<cy[S[p].y1];i++){
					int k=vy[S[p].y1][i];
					if(k==p)continue;
					if(intersect(nr,S[k])){
						T[c]=k;
						c++;
					}
				}
				double d;
				d=calc_score(nr,R[p])-calc_score(S[p],R[p]);
				for(int i=0;i<c;i++){
					rect tmprec=S[T[i]];
					tmprec.y2--;
					if(tmprec.y2<=Y[T[i]]){
                      	d-=334334334;
                    }
					else{
						d+=calc_score(tmprec,R[T[i]])-calc_score(S[T[i]],R[T[i]]);
					}
				}
				if(c==0&&d<0&&S[p].y2!=Y[p]+1){
					nr.y2--;
					while(true){
						if(nr.y1==0||nr.y2==Y[p]+1)break;
						bool f=true;
						nr.y1--;
						nr.y2--;
						for(int i=0;i<cy[nr.y1+1];i++){
							if(vy[nr.y1+1][i]==p)continue;
							if(intersect(nr,S[vy[nr.y1+1][i]])){
								f=false;
                              	break;
							}
						}
						if(!f){
							nr.y1++;
							nr.y2++;
							break;
						}
					}
					rect_update(S[p],nr,p);
				}
				else if(d>=-(double)_/C1){
					score+=d;
					rect_update(S[p],nr,p);
					for(int i=0;i<c;i++){
						rect tmprec=S[T[i]];
						tmprec.y2--;
						rect_update(S[T[i]],tmprec,T[i]);
					}
				}
			}
			if(score>maxscore){
				maxscore=score;
				for(int i=0;i<N;i++){
					resA[i]=S[i].x1;
					resB[i]=S[i].y1;
					resC[i]=S[i].x2;
					resD[i]=S[i].y2;
				}
			}
		}
		if(_%4==2){
			p=xor128()%N;
			//+x
			if(S[p].x2!=10000){
				rect nr=S[p];
				nr.x2++;
				if(_%M==0&&S[p].y1!=Y[p]){
					nr.y1++;
				}
				if(_%M==1&&S[p].y2!=Y[p]+1){
					nr.y2--;
				}
				int c=0;
				for(int i=0;i<cx[S[p].x2];i++){
					int k=vx[S[p].x2][i];
					if(k==p)continue;
					if(intersect(nr,S[k])){
						T[c]=k;
						c++;
					}
				}
				double d;
				d=calc_score(nr,R[p])-calc_score(S[p],R[p]);
				for(int i=0;i<c;i++){
					rect tmprec=S[T[i]];
					tmprec.x1++;
					if(tmprec.x1>X[T[i]]){
                      	d-=334334334;
                    }
					else{
						d+=calc_score(tmprec,R[T[i]])-calc_score(S[T[i]],R[T[i]]);
					}
				}
				if(c==0&&d<0&&S[p].x1!=X[p]){
					nr.x1++;
					while(true){
						if(nr.x1==X[p]||nr.x2==10000)break;
						bool f=true;
						nr.x1++;
						nr.x2++;
						for(int i=0;i<cx[nr.x2-1];i++){
							if(vx[nr.x2-1][i]==p)continue;
							if(intersect(nr,S[vx[nr.x2-1][i]])){
								f=false;
                              	break;
							}
						}
						if(!f){
							nr.x1--;
							nr.x2--;
							break;
						}
					}
					rect_update(S[p],nr,p);
				}
				else if(d>=-(double)_/C1){
					score+=d;
					rect_update(S[p],nr,p);
					for(int i=0;i<c;i++){
						rect tmprec=S[T[i]];
						tmprec.x1++;
						rect_update(S[T[i]],tmprec,T[i]);
					}
				}
			}
			if(score>maxscore){
				maxscore=score;
				for(int i=0;i<N;i++){
					resA[i]=S[i].x1;
					resB[i]=S[i].y1;
					resC[i]=S[i].x2;
					resD[i]=S[i].y2;
				}
			}
		}
		if(_%4==3){
			p=xor128()%N;
			//+y
			if(S[p].y2!=10000){
				rect nr=S[p];
				nr.y2++;
				if(_%M==0&&S[p].x1!=X[p]){
					nr.x1++;
				}
				if(_%M==1&&S[p].x2!=X[p]+1){
					nr.x2--;
				}
				int c=0;
				for(int i=0;i<cy[S[p].y2];i++){
					int k=vy[S[p].y2][i];
					if(k==p)continue;
					if(intersect(nr,S[k])){
						T[c]=k;
						c++;
					}
				}
				double d;
				d=calc_score(nr,R[p])-calc_score(S[p],R[p]);
				for(int i=0;i<c;i++){
					rect tmprec=S[T[i]];
					tmprec.y1++;
					if(tmprec.y1>Y[T[i]]){
                      	d-=334334334;
                    }
					else{
						d+=calc_score(tmprec,R[T[i]])-calc_score(S[T[i]],R[T[i]]);
					}
				}
				if(c==0&&d<0&&S[p].y1!=Y[p]){
					nr.y1++;
					while(true){
						if(nr.y1==Y[p]||nr.y2==10000)break;
						bool f=true;
						nr.y1++;
						nr.y2++;
						for(int i=0;i<cy[nr.y2-1];i++){
							if(vy[nr.y2-1][i]==p)continue;
							if(intersect(nr,S[vy[nr.y2-1][i]])){
								f=false;
                              	break;
							}
						}
						if(!f){
							nr.y1--;
							nr.y2--;
							break;
						}
					}
					rect_update(S[p],nr,p);
				}
				else if(d>=-(double)_/C1){
					score+=d;
					rect_update(S[p],nr,p);
					for(int i=0;i<c;i++){
						rect tmprec=S[T[i]];
						tmprec.y1++;
						rect_update(S[T[i]],tmprec,T[i]);
					}
				}
			}
			if(score>maxscore){
				maxscore=score;
				for(int i=0;i<N;i++){
					resA[i]=S[i].x1;
					resB[i]=S[i].y1;
					resC[i]=S[i].x2;
					resD[i]=S[i].y2;
				}
			}
		}
	}
	for(int i=0;i<N;i++){
		printf("%d %d %d %d\n",resA[i],resB[i],resC[i],resD[i]);
	}
}

 
