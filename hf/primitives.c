#include "primitives.h"
#include <stdio.h>
#include "basis.h"
prim build_primitive(float coef,float alpha,int l, int m, int n){
    float norm = sqrt(pow(2 * alpha/ M_PI,1.5f) * powf(4.0f * alpha,l+m+n)/
                    (double_factorial(2*l-1) * double_factorial(2* m-1) * double_factorial(2*n-1)));

    prim p = {
        norm,coef,alpha,l,m,n
    };
    return p;
}
int build_cgto(int at_num, atoms  * ats, cgto * cgtos, prim * primitives){
    //number of cgtos must be precalculated before passing array or crash because no malloc
    //returns number of cgtos
    prim * p = primitives;
    int orb = 0;

    for (int a = 0; a < at_num; a++){
            //printf("%i",(int)p);
        int num_orbs = ORBS[ats[a].type-1];
        p[0] = build_primitive(D1S[0],AS1[ats[a].type - 1][0],0,0,0);
        p[1] = build_primitive(D1S[1],AS1[ats[a].type - 1][1],0,0,0);
        p[2] = build_primitive(D1S[2],AS1[ats[a].type - 1][2],0,0,0);
        cgtos[orb] = (cgto){{ats[a].pos[0],ats[a].pos[1],ats[a].pos[2]},p,3};
        p+=3;
        orb+=1;

        if (num_orbs == 5){
            p[0] = build_primitive(D2S[0],AS2[ats[a].type - 1][0],0,0,0);
            p[1] = build_primitive(D2S[1],AS2[ats[a].type - 1][1],0,0,0);
            p[2] = build_primitive(D2S[2],AS2[ats[a].type - 1][2],0,0,0);
            cgtos[orb] = (cgto){{ats[a].pos[0],ats[a].pos[1],ats[a].pos[2]},p,3};
            p+=3;
            orb+=1;
            p[0] = build_primitive(D2P[0],AS2[ats[a].type - 1][0],1,0,0);
            p[1] = build_primitive(D2P[1],AS2[ats[a].type - 1][1],1,0,0);
            p[2] = build_primitive(D2P[2],AS2[ats[a].type - 1][2],1,0,0);
            cgtos[orb] = (cgto) {{ats[a].pos[0],ats[a].pos[1],ats[a].pos[2]},p,3};
            p+=3;
            orb+=1;
            p[0] = build_primitive(D2P[0],AS2[ats[a].type - 1][0],0,1,0),
            p[1] = build_primitive(D2P[1],AS2[ats[a].type - 1][1],0,1,0),
            p[2] = build_primitive(D2P[2],AS2[ats[a].type - 1][2],0,1,0),
            cgtos[orb]  = (cgto){{ats[a].pos[0],ats[a].pos[1],ats[a].pos[2]},p,3};
            p+=3;
            orb+=1;
            p[0] = build_primitive(D2P[0],AS2[ats[a].type - 1][0],0,0,1),
            p[1] = build_primitive(D2P[1],AS2[ats[a].type - 1][1],0,0,1),
            p[2] = build_primitive(D2P[2],AS2[ats[a].type - 1][2],0,0,1),
            cgtos[orb] = (cgto){{ats[a].pos[0],ats[a].pos[1],ats[a].pos[2]},p,3};
            p+=3;
            orb+=1;
        }
    }
    return orb;
}