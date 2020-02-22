#include <stdio.h>
#include <stdlib.h>
#include <math.h>


typedef struct
{
    unsigned char r,g,b;
} pixel;

typedef struct
{
    unsigned int W, H, dim_img, dim_padding;
    unsigned char header[54];
    pixel *vector;

} imagine;

typedef struct
{
    unsigned int x,y,W,H;
    double indice_corelare;
    unsigned sablon_folosit;
} fereastra;

typedef struct
{
    fereastra *v_ferestre;
    unsigned int numar_ferestre;

} vector_cu_ferestre;

typedef struct
{
    unsigned int W, H, dim_img, dim_padding, numar_pixeli;
    unsigned char header[54];
    pixel* vector;
    double deviatie_standard;
    double medie;

} info_sablon;


imagine citire_imagine(char * nume_fisier)
{
    FILE * fin = fopen(nume_fisier, "rb");
    if (fin == NULL)
    {
        printf("eroare");
        exit(1);
    }

    imagine img;
    fread(img.header, sizeof(unsigned char), 54, fin);

    fseek(fin, 2, SEEK_SET);
    fread(&img.dim_img, sizeof(img.dim_img), 1, fin);
    //printf("%u", img.dim_img);
    fseek(fin, 18, SEEK_SET);
    fread(&img.W, sizeof(img.W), 1, fin);
    //printf("\n %u", img.W);
    fseek(fin, 22, SEEK_SET);
    fread(&img.H, sizeof(img.H), 1, fin);
    //printf("\n %u",img.H);
    if (img.W % 4 != 0)
        img.dim_padding = 4 - 3 * img.W % 4;
    else
        img.dim_padding = 0;
    //  printf ("\n %u", img.dim_padding);

    int i, j, k;
    unsigned char c;
    fseek(fin, 54, SEEK_SET);
    img.vector = (pixel *)calloc(img.W * img.H, sizeof(pixel));

    for(i = img.H - 1; i >=0; i--)
    {
        for(j = 0; j < img.W; j++)
        {
            fread(&(img.vector[img.W * i + j].b), 1, 1, fin);
            fread(&(img.vector[img.W * i + j].g), 1, 1, fin);
            fread(&(img.vector[img.W * i + j].r), 1, 1, fin);
        }

        for(k = 0; k < img.dim_padding; k++)
            fread(&c, 1, 1, fin);
    }

    fclose(fin);

    return img;
}
    unsigned int xorshift ( unsigned int numar_generat_anterior)
{
    unsigned int numar_curent = numar_generat_anterior;
    numar_curent = numar_curent ^ numar_curent << 13;
    numar_curent = numar_curent ^ numar_curent >> 17;
    numar_curent = numar_curent ^ numar_curent << 5;
    return numar_curent;

}

pixel pixel_xor_constanta (pixel pixel1, unsigned int x)
{
    pixel rezultat;
    rezultat = pixel1;
    unsigned char x0,x1,x2;
    x0= x & 255;
    x=x>>8;
    x1=x&255;
    x=x>>8;
    x2=x&255;
    rezultat.b= x0^rezultat.b;
    rezultat.g=x1^rezultat.g;
    rezultat.r=x2^rezultat.r;
    return rezultat;

}
pixel pixeli_xor_pixeli (pixel pixel1, pixel pixel2)
{
    pixel rezultat;
    rezultat.b= pixel1.b^pixel2.b;
    rezultat.g= pixel1.g^pixel2.g;
    rezultat.r=pixel1.r^pixel2.r;
    return rezultat;
}
pixel *  criptare (imagine img, char * fisier_cheie)
{
    unsigned int a, SV;

    FILE *fin;
    fin = fopen (fisier_cheie, "r");
    fscanf(fin, "%u %u", &a, &SV);




    unsigned  int *R, p, i;
    R = (unsigned int *) malloc (sizeof(unsigned int) * 2 * img.H * img.W );
    p = a;
    R[0] = a;
    for(i = 1; i <= 2 * img.H * img.W - 1; i++)
    {
        p  = xorshift(p);
        R[i] = p;

    }

    unsigned int *permutare;
    permutare  = (unsigned int *)calloc(img.H * img.W, sizeof(unsigned int));
    int k = img.H * img.W - 1;
    int r;

    for (i = 1; i <= k; i++)
        permutare[i] = i;

    int j, aux;
    int s = 1;
    for (j = k; j >= 0; j--)
    {
        r = R[s] % (j + 1);
        s++;

        aux = permutare[r];
        permutare[r] = permutare[j];
        permutare[j] = aux;
    }


    imagine img2;
    img2 = img;

    img2.vector = (pixel *)malloc(sizeof(pixel) * img.H * img.W);

    for (i = 0; i <= k; i++)
    {

        img2.vector[permutare[i]] = img.vector[i];

    }
    pixel *C;
    C = (pixel *)malloc(sizeof(pixel) * img.W * img.H);


    for (i = 0; i <= k ; i++)
    {
        if (i == 0)
            C[0] = pixel_xor_constanta(pixel_xor_constanta(img2.vector[i], SV), R[img2.W * img2.H]);


        else

            C[i]=  pixel_xor_constanta(pixeli_xor_pixeli(C[i-1], img2.vector[i]), R[img2.H * img2.W + i]);


    }
    fclose(fin);

    return C;
}
imagine decriptare (imagine img, char * fisier_cheie)
{
    unsigned int a, SV;

    FILE *fin;
    fin = fopen (fisier_cheie, "r");
    fscanf(fin, "%u %u", &a, &SV);


    unsigned  int *R, p, i;
    R = (unsigned int *) malloc (sizeof(unsigned int) * 2 * img.H * img.W );
    p = a;
    R[0] = a;
    for(i = 1; i <= 2 * img.H * img.W - 1; i++)
    {
        p  = xorshift(p);
        R[i] = p;

    }

    unsigned int *permutare;
    permutare  = (unsigned int *)calloc(img.H * img.W, sizeof(unsigned int));
    int k = img.H * img.W - 1;
    int r;
    for (i = 1; i <= k; i++)
        permutare[i] = i;

    int j, aux;
    int s = 1;
    for (j = k; j >= 0; j--)
    {
        r = R[s] % (j + 1);
        s++;

        aux = permutare[r];
        permutare[r] = permutare[j];
        permutare[j] = aux;
    }

    unsigned int *permutareinv;
    permutareinv  = (unsigned int *)calloc(img.H * img.W, sizeof(unsigned int));
    for(j=k; j>=0; j--)
    {
        permutareinv[permutare[j]] = j;
    }


    imagine Cprim;
    Cprim = img;
    Cprim.vector = (pixel *)malloc(sizeof(pixel) * img.H * img.W);

    for (i = 0; i <=k  ; i++)

        if (i == 0)

            Cprim.vector[i] = pixel_xor_constanta(pixel_xor_constanta(img.vector[0], SV), R[img.W * img.H]);

        else
            Cprim.vector[i]=  pixel_xor_constanta(pixeli_xor_pixeli(img.vector[i-1], img.vector[i]), R[img.H * img.W + i]);

    imagine D;
    D = Cprim;
    D.vector =(pixel *)malloc(sizeof(pixel) * img.H * img.W);

    for (i = 0 ; i <= k ; i++)
    {
        D.vector[permutareinv[i]] = Cprim.vector[i];
    }

    fclose(fin);

    return D;
}
void chi_patrat (imagine img, char *nume_fisier)
{
    pixel *f;
    f = (pixel *) malloc (sizeof(pixel)* img.H*img.W);
    FILE *fin;
    int i;
    fin = fopen(nume_fisier,"rb");
    unsigned int fblue[256], fgreen[256], fred[256];


   for(i=0 ; i <=255 ; i++)
   {
       fblue[i] = 0;
       fred[i]= 0;
       fgreen[i]=0;
   }
      fseek(fin, 54, SEEK_SET);
    for(i=0 ; i < img.H * img. W ; i++)
    {
        fread(&(img.vector[i].b),1,1,fin);
        fblue[img.vector[i].b]++;
        fread(&(img.vector[i].g),1,1,fin);
        fgreen[img.vector[i].g]++;
        fread(&(img.vector[i].r),1,1,fin);
        fred[img.vector[i].r]++;

    }

    double chi_b,chi_r,chi_g;
    chi_b=0;
    chi_r=0;
    chi_g=0;
    double frecventa_teoretica;
    frecventa_teoretica = (1.0 * (img.H * img.W))/256.0;
    for(i=0; i <=255 ; i++)
    {
        chi_b = chi_b + ((fblue[i] - frecventa_teoretica)* (fblue[i] - frecventa_teoretica))/frecventa_teoretica;
        chi_g= chi_g + ((fgreen[i] - frecventa_teoretica)* (fgreen[i] - frecventa_teoretica))/frecventa_teoretica;
        chi_r = chi_r + ((fred[i] - frecventa_teoretica)* (fred[i] - frecventa_teoretica))/frecventa_teoretica;

    }

    printf("\n  B %.2f ",chi_r);
    printf("\n  G %.2f",chi_g);
    printf("\n  R %.2f",chi_b);
}

void afisareImagine(char * nume_fisier, imagine img)
{

    FILE * fout = fopen(nume_fisier, "wb");
    int i, j, k;
    unsigned char c = 0;
    for (i = 0; i < 54; ++i)
    {
        fwrite(&(img.header[i]), sizeof(unsigned char), 1, fout);
    }

    for (i = img.H - 1; i >= 0; --i)
    {
        for (j = 0; j < img.W; ++j)
        {
            fwrite(&(img.vector[img.W * i + j].b), 1, 1, fout);
            fwrite(&(img.vector[img.W * i + j].g), 1, 1, fout);
            fwrite(&(img.vector[img.W * i + j].r), 1, 1, fout);
        }

        for (k = 0; k < img.dim_padding; ++k)
            fwrite(&c, 1, 1, fout);
    }

    fclose(fout);
}


//fin

imagine grayscale ( imagine img)
{
    int i;
    char aux;
    for(i = 0; i < img.H * img.W ; i++)

    {
        aux = 0.299*img.vector[i].r + 0.587*img.vector[i].g + 0.114*img.vector[i].b;
        img.vector[i].r = img.vector[i].g = img.vector[i].b = aux;
    }
    return img;
}

info_sablon grayscale_sablon ( info_sablon img)
{
    int i;
    char aux;
    for(i = 0; i < img.H * img.W ; i++)

    {
        aux = 0.299*img.vector[i].r + 0.587*img.vector[i].g + 0.114*img.vector[i].b;
        img.vector[i].r = img.vector[i].g = img.vector[i].b = aux;
    }
    return img;
}
info_sablon citire_sablon (char *nume_sablon)
{
    FILE * fin = fopen(nume_sablon, "rb");
    if (fin == NULL)
    {
        printf("eroare");
        exit(1);
    }

    info_sablon img;
    fread(img.header, sizeof(unsigned char), 54, fin);

    fseek(fin, 2, SEEK_SET);
    fread(&img.dim_img, sizeof(img.dim_img), 1, fin);
    //printf("%u", img.dim_img);
    fseek(fin, 18, SEEK_SET);
    fread(&img.W, sizeof(img.W), 1, fin);
    //printf("\n %u", img.W);
    fseek(fin, 22, SEEK_SET);
    fread(&img.H, sizeof(img.H), 1, fin);
    //printf("\n %u",img.H);
    if (img.W % 4 != 0)
        img.dim_padding = 4 - 3 * img.W % 4;
    else
        img.dim_padding = 0;
    //  printf ("\n %u", img.dim_padding);

    int i, j, k;
    unsigned char c;
    fseek(fin, 54, SEEK_SET);
    img.vector = (pixel *)calloc(img.W * img.H, sizeof(pixel));
    img.numar_pixeli = img.H * img.W;
    double Suma = 0.0f, suma_deviatie = 0.0f;

    for(i = img.H - 1; i >=0; i--)
    {
        for(j = 0; j < img.W; j++)
        {
            fread(&(img.vector[img.W * i + j].b), 1, 1, fin);
            fread(&(img.vector[img.W * i + j].g), 1, 1, fin);
            fread(&(img.vector[img.W * i + j].r), 1, 1, fin);
            Suma = Suma + img.vector[i* img.W + j].r;
        }

        for(k = 0; k < img.dim_padding; k++)
            fread(&c, 1, 1, fin);
    }

    img.medie = Suma / (img.H * img.W);

    for(i = img.H - 1; i >=0; i--)
        for(j = 0; j < img.W; j++)
            suma_deviatie = suma_deviatie + (img.vector[i* img.W + j].r - img.medie) * (img.vector[i* img.W + j].r - img.medie);

    img.deviatie_standard = sqrt(suma_deviatie/(img.numar_pixeli-1));


    fclose(fin);

    return img;
}

double medie_fereastra (imagine frst, info_sablon sablon)
{
    double media = 0.0f, suma = 0.0f;
    int i,j;

    for(i= sablon.H -1 ; i >= 0 ; i --)
        for(j=0 ; j< sablon.W; j++)
            suma = suma + frst.vector[i * sablon.W + j].r;
    media = suma / (sablon.H  * sablon.W);

    return media;
}
double deviatie_fereastra ( imagine img, info_sablon sablon)
{
    int i,j;
    double media = 0.0f, suma = 0.0f, deviatie = 0.0f;
    media = medie_fereastra(img,sablon);
    for(i= sablon.H -1 ; i >= 0 ; i --)
        for(j=0 ; j< sablon.W; j++)
            suma = suma + (img.vector[ i * sablon.W +j].r - media)*(img.vector[i * sablon.W + j].r - media);


    deviatie = sqrt(suma/(sablon.H * sablon.W - 1));
    return deviatie;
}

double corr(info_sablon sablon, imagine frst)
{
    int i,j;
    double suma = 0.0f, corelatie = 0.0f;
    double deviat_fereastra= deviatie_fereastra(frst,sablon);
    double med_fereastra = medie_fereastra(frst,sablon);

    for(i = sablon.H -1 ; i >=0 ; i--)
    {
        for(j=0; j < sablon.W ; j++)
        {
            suma = suma + ((frst.vector[i*sablon.W + j].r - med_fereastra) * (sablon.vector[ sablon.W * i + j].r - sablon.medie))/(deviat_fereastra*sablon.deviatie_standard);

        }

    }
    corelatie = suma/(sablon.H * sablon.W);
    return corelatie;
}
vector_cu_ferestre repartizare(imagine img_principala, info_sablon sablon)
{
    double prag = 0.50;
    double p;
    vector_cu_ferestre vector;
    vector.v_ferestre = (fereastra *) malloc (sizeof(fereastra) * img_principala.H * img_principala.W);
    int k=0;
    imagine fr;
    fr.vector = (pixel *)calloc(sablon.W * sablon.H, sizeof(pixel));
    int i,j, linie, coloana;
    for (i = 0; i <= img_principala.H - sablon.H; i++)
    {
        for (j = 0; j <= img_principala.W - sablon.W; j++)
        {
            for (linie = i; linie < i + sablon.H; linie++)
            {
                for (coloana = j; coloana < j + sablon.W; coloana++)
                {
                    fr.vector[(linie - i)* sablon.W + coloana - j ] = img_principala.vector[linie* img_principala.W + coloana];
                }
            }

            p = corr(sablon,fr);

            if(p >= prag)
            {
                k++;
                vector.v_ferestre[k].H = sablon.H;
                vector.v_ferestre[k].W = sablon.W;
                vector.v_ferestre[k].indice_corelare = p;
                vector.v_ferestre[k].x = j;
                vector.v_ferestre[k].y = i;
            }
        }
    }

    vector.numar_ferestre = k;
    return vector;
}

pixel get_culoare(int k)
{
    switch(k)
    {
    case 0:
        return (pixel)
        {
            255, 0, 0
        };
    case 1:
        return (pixel)
        {
            255, 255, 0
        };
    case 2:
        return (pixel)
        {
            0, 255, 0
        };
    case 3:
        return (pixel)
        {
            0, 255, 255
        };
    case 4:
        return (pixel)
        {
            255, 0, 255
        };
    case 5:
        return (pixel)
        {
            0, 0, 255
        };
    case 6:
        return (pixel)
        {
            192, 192, 192
        };
    case 7:
        return (pixel)
        {
            255, 140, 0
        };
    case 8:
        return (pixel)
        {
            128, 0, 128
        };
    case 9:
        return (pixel)
        {
            128, 0, 0
        };
    default:
        return (pixel)
        {
            255, 255, 255
        };
    }
}

imagine colorare_img (imagine imagine_principala, info_sablon sablon,  vector_cu_ferestre *vector)
{
    int i,j;

    pixel culoare;
    for (i = 0 ; i < vector->numar_ferestre ; i ++)
    {
        culoare=get_culoare(vector->v_ferestre[i].sablon_folosit);

        for (j = vector->v_ferestre[i].x ; j <= vector->v_ferestre[i].x +  sablon.W ; j++)
        {
            imagine_principala.vector[vector->v_ferestre[i].y * imagine_principala.W + j]=culoare;
            imagine_principala.vector[(vector->v_ferestre[i].y + sablon.H) * imagine_principala.W + j]=culoare;
        }

        for (j = vector->v_ferestre[i].y ; j <= vector->v_ferestre[i].y +  vector->v_ferestre[i].H; j++)
        {
            imagine_principala.vector[j * imagine_principala.W + vector->v_ferestre[i].x] =(pixel)culoare;
            imagine_principala.vector[j * imagine_principala.W + vector->v_ferestre[i].x + sablon.W] =(pixel)culoare;
        }

    }

    return imagine_principala;
}

int cmp_detectii (const void *a, const void *b)
{
    fereastra *fa = (fereastra *) a;
    fereastra *fb = (fereastra *) b;

    if(fa->indice_corelare < fb ->indice_corelare)
        return 1;
    return -1;

}

double max (double a, double b)
{
    return (a > b) ? a : b;
}

double min (double a, double b)
{
    return (a < b) ? a : b;
}
int suprapunere (fereastra f1, fereastra f2)
{
    double suprafata1, suprafata2, suprafata3, suprafata4,ratia;
    suprafata1 = f1.W * f1.H * 1.0;

    suprafata2 = f2.H * f2.W * 1.0;

    suprafata3 = max(0, min(f2.x +f2.W, f1.x +f1.W) - max(f2.x, f1.x)) * max(0, min(f2.y +f2.H, f1.y +f2.H) - max(f2.y, f1.y)) * 1.0;
    suprafata4 = suprafata1 + suprafata2 - suprafata3;
    ratia = suprafata3 / suprafata4;

    if(ratia > 0.2)
        return 1;
    return 0;
}


void eliminare_suprapuneri ( vector_cu_ferestre *vector )
{
    int i, j,numar, l;
    qsort(vector->v_ferestre, vector->numar_ferestre,sizeof(fereastra),cmp_detectii);

    numar = vector->numar_ferestre;
    for(i=0 ; i < numar; i++)
    {
        for(j = i +1 ; j < numar ; j++ )
        {
            if(suprapunere(vector->v_ferestre[i],vector->v_ferestre[j]) == 1 )
            {
                for( l = j ; l < numar ; l++ )
                    vector->v_ferestre[l] = vector->v_ferestre[l+1];
                numar--;
            }
        }
    }

    vector->numar_ferestre = numar;

}
int main()
{
    imagine img = citire_imagine("peppers.bmp");

    pixel * C = criptare(img, "key.txt");
    imagine crip = img;
    crip.vector = C;
    afisareImagine("criptare.bmp", crip);

    imagine decriptata = decriptare(crip, "key.txt");
    afisareImagine("decriptare.bmp",decriptata);
    printf("\n Chi-squared test on RGB channels for peppers.bmp:");
    chi_patrat(img,"P.bmp");
    printf("\n Chi-squared test on RGB channels for criptare.bmp:");
    chi_patrat(crip,"criptare.bmp");


    int i;
    imagine img_principala = citire_imagine("test.bmp");
    vector_cu_ferestre vector, vector_mare;
    imagine grayscale_test = grayscale(img_principala);
    afisareImagine("grayscale_test.bmp",grayscale_test);

    imagine colorare = img_principala;
    pixel culoare;
    info_sablon sablon;

    char* sablon_cale = (char*)malloc(128);

    vector_mare.numar_ferestre = 0;
    vector_mare.v_ferestre = (fereastra*)malloc(sizeof(fereastra) * 100000);
    info_sablon sab;

    for(unsigned i = 0; i < 10; ++i)
    {
        sprintf(sablon_cale, "cifra%u.bmp", i);

        sablon = citire_sablon(sablon_cale);
       sab = grayscale_sablon(sablon);
        vector = repartizare(colorare, sab);

        if(vector_mare.numar_ferestre + vector.numar_ferestre > 100000)
        {
            printf("prea multe ferestre\n");
            return 0;
        }

        for(unsigned j = vector_mare.numar_ferestre; j < vector_mare.numar_ferestre + vector.numar_ferestre; ++j)
        {
            vector_mare.v_ferestre[j] = vector.v_ferestre[j - vector_mare.numar_ferestre];
            vector_mare.v_ferestre[j].sablon_folosit = i;
        }
        vector_mare.numar_ferestre += vector.numar_ferestre;
        free(vector.v_ferestre);
    }


    eliminare_suprapuneri(&vector_mare);
    colorare = colorare_img(colorare, sab, &vector_mare);

    afisareImagine("elim.bmp",colorare);
    free(vector_mare.v_ferestre);
    free(sablon_cale);
    return 0;
}


