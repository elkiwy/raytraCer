#include "perlin.h"




void permute(int* p, int n){
    for(int i=0;i<n;++i){
        int target = random_int(0, i);
        int tmp = p[i];
        p[i] = p[target];
        p[target] = tmp;
    }
}


int* perlin_generate_perm(){
    int* p = malloc(sizeof(int)*POINT_COUNT);
    for(int i=0;i<POINT_COUNT;++i){p[i] = i;}
    permute(p, POINT_COUNT);
    return p;
}

double trilinear_interp(double c[2][2][2], double u, double v, double w){
    double accum = 0.0;
    for (int i=0;i<2;i++){
        for (int j=0;j<2;j++){
            for (int k=0;k<2;k++){
                accum += (i*u + (1-i)*(1-u)) * (j*v + (1-j)*(1-v)) * (k*w + (1-k)*(1-w)) * c[i][j][k];
            }
        }
    }
    return accum;
}


double perlin_interp(vec3 c[2][2][2], double u, double v, double w){
    double uu = u*u*(3-2*u);// Hermitan smoothing
    double vv = v*v*(3-2*v);
    double ww = w*w*(3-2*w);
    double accum = 0.0;
    for (int i=0;i<2;i++){
        for (int j=0;j<2;j++){
            for (int k=0;k<2;k++){
                vec3 weight_v = {u-i, v-j, w-k};
                accum += (i*uu + (1-i)*(1-uu))
                       * (j*vv + (1-j)*(1-vv))
                       * (k*ww + (1-k)*(1-ww))
                       * vec3c_dot(c[i][j][k], weight_v);
            }
        }
    }
    return accum;
}











perlin* perlin_init(){
    perlin* p = malloc(sizeof(perlin));
    p->ranvec = malloc(sizeof(vec3)*POINT_COUNT);
    for (int i=0;i<POINT_COUNT; ++i){
        p->ranvec[i] = vec3c_unit(vec3_random_scaled(-1, 1));
    }
    p->perm_x = perlin_generate_perm();
    p->perm_y = perlin_generate_perm();
    p->perm_z = perlin_generate_perm();
    return p;
}




void perlin_free(perlin* p){
    free(p->ranvec);
    free(p->perm_x);
    free(p->perm_y);
    free(p->perm_z);
    free(p);
}





double perlin_noise(perlin* p, point3 point){
    double u = point.x - floor(point.x);
    double v = point.y - floor(point.y);
    double w = point.z - floor(point.z);
    int i = (int)floor(point.x);
    int j = (int)floor(point.y);
    int k = (int)floor(point.z);

    vec3 c[2][2][2];

    for (int di=0; di < 2; di++){
        for (int dj=0; dj < 2; dj++){
            for (int dk=0; dk < 2; dk++){
                c[di][dj][dk] = p->ranvec[
                    p->perm_x[(i+di) & 255] ^
                    p->perm_y[(j+dj) & 255] ^
                    p->perm_z[(k+dk) & 255]
                ];
            }
        }
    }

    return perlin_interp(c, u, v, w);
}




double turb_rec(perlin* perl, point3 p, int depth){
    double accum = 0.0;
    point3 temp_p = p;
    double weight = 1.0;
    for (int i=0; i<depth; i++){
        accum += weight * perlin_noise(perl, temp_p);
        weight *= 0.5;
        temp_p = vec3c_mul_k(temp_p, 2);
    }
    return fabs(accum);
}
double perlin_turb(perlin* perl, point3 p){
    return turb_rec(perl, p, 7);
}
