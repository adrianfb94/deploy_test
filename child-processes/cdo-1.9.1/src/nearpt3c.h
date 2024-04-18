//typedef unsigned short int Coord_T;
//typedef short int Coord_T;
//typedef unsigned int Coord_T;
typedef float Coord_T;

#define NPT3SFACT 32000
#define NPT3SCALE(x) (0.5+(x+1)*NPT3SFACT)

#if defined(__cplusplus)
extern "C" {
#endif
void *nearpt3_preprocess(const int nfixpts, Coord_T **pts);
int nearpt3_query(void *g, const Coord_T *q);
void nearpt3_destroy(void *g);
#if defined(__cplusplus)
}
#endif
