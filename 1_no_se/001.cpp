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

int N,X[200],Y[200],R[200],ord[200];
rect S[200];

inline bool include(rect r,int p,int q){
	return r.x1<=p&&p<=r.x2&&r.y1<=q&&q<=r.y2;
}

inline bool across(rect r1,rect r2){
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

int main(){
	scanf("%d",&N);
	for(int i=0;i<N;i++){
		scanf("%d%d%d",&X[i],&Y[i],&R[i]);
	}
	for(int i=0;i<N;i++){
		S[i].x1=X[i];
		S[i].y1=Y[i];
		S[i].x2=X[i]+1;
		S[i].y2=Y[i]+1;
		ord[i]=i;
	}
	for(int _=100000;_>=0;_--){
		int p,q;
		p=xor128()%10000;
		q=xor128()%10000;
		random_shuffle(ord,ord+N);
		for(int i=0;i<N;i++){
			int t=ord[i];
			rect nr=newrect(S[t],X[t],Y[t],p,q);
			if(nr.area()>R[t]*2)continue;
			if(-(abs(1-(double)min(R[t],nr.area())/max(R[t],nr.area()))-abs(1-(double)min(R[t],S[t].area())/max(R[t],S[t].area())))>=-(double)_/1000000){
				bool f=true;
				for(int j=0;j<N;j++){
					if(j==t)continue;
					if(across(nr,S[j])){
						f=false;
						break;
					}
				}
				if(f){
					S[t]=nr;
					break;
				}
			}
		}
	}
	for(int i=0;i<N;i++){
		printf("%d %d %d %d\n",S[i].x1,S[i].y1,S[i].x2,S[i].y2);
	}
}
