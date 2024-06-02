#include "lodepng.c"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define pi 2*acos(0)

typedef struct pixel
{
	int R;
	int G;
	int B;
	int alpha;
} pixel;

struct DUS;

typedef struct Elem
{
    struct DUS* head;
    pixel* data;
    struct Elem* next;
} Elem;

typedef struct DUS
{
    int size;
    int color_status;
    struct Elem* beg, *end;
} DUS;


char* load_png_file  (const char *filename, int *width, int *height) {
	unsigned char *image = NULL;
	int error = lodepng_decode32_file(&image, width, height, filename);
	if (error) {
		printf("error %u: %s\n", error, lodepng_error_text(error));
		return NULL;
	}
	return (image);
}

int code_png_file  (const char *filename, unsigned char *image, int width, int height) {
	int error = lodepng_encode32_file(filename, image, width, height);
	if (error) {
		printf("error %u: %s\n", error, lodepng_error_text(error));
	}
	return (error);
}

pixel** create_image(int h, int w, char* picture)
{
    int i, j, k=0;
    pixel** Pict = malloc(h*sizeof(pixel*));
    for (i=0; i<h; i++) Pict[i] = malloc(w*sizeof(pixel));
    if (picture==NULL) for (i=0; i<h; i++) for (j=0; j<w; j++) {
        Pict[i][j].R = 0;
        Pict[i][j].G = 0;
        Pict[i][j].B = 0;
        Pict[i][j].alpha = 0;
        k += 4;
    }
    else for (i=0; i<h; i++) for (j=0; j<w; j++) {
        Pict[i][j].R = picture[k+0];
        Pict[i][j].G = picture[k+1];
        Pict[i][j].B = picture[k+2];
        Pict[i][j].alpha = picture[k+3];
        k += 4;
    }
    return Pict;
}

void clear_image(int h, int w, pixel** Pict)
{
    int i=0;
    for (i; i<h; i++) free(Pict[i]);
    free(Pict);
    return;
}

pixel** get_image(char* filename, int* height, int* width)
{
    int h=0, w=0;
    char* picture = load_png_file(filename, &w, &h);
	if (picture == NULL) {
		printf("Error by reading the initial file.\n");
		return NULL;
	}
	pixel** Pict = create_image(h, w, picture);
	*height = h;
	*width = w;
	return Pict;
}

void put_image(int h, int w, pixel** Pict, char* filename)
{
    int i, j;
    char* picture = malloc(4*h*w*sizeof(char));
    for (i=0; i<h; i++) for (j=0; j<w; j++) {
        picture[4*(w*i+j)+0] = Pict[i][j].R;
        picture[4*(w*i+j)+1] = Pict[i][j].G;
        picture[4*(w*i+j)+2] = Pict[i][j].B;
        picture[4*(w*i+j)+3] = Pict[i][j].alpha;
    }
    int err = code_png_file(filename, picture, w, h);
	if (err) printf("Error by writing in the outer file.\n");
	free(picture);
    return;
}

void image_to_BW(int h, int w, pixel** Pict, int border_value)
{
    int i, j, res;
    for (i=0; i<h; i++) for (j=0; j<w; j++) {
        res = (Pict[i][j].R + Pict[i][j].G + Pict[i][j].B) / 3;
        if (res < 0) res += 255;
        if (res < border_value) res = 0;
        Pict[i][j].R = res;
        Pict[i][j].G = res;
        Pict[i][j].B = res;
        Pict[i][j].alpha = 255;
    }
    return;
}

void Gauss_operator(int h, int w, pixel** Pict, pixel** Pict_new, double sigma)
{
    int i, j, p, q;
    double res, G[5][5];
    for (i=-2; i<3; i++) for (j=-2; j<3; j++) {
        G[i+2][j+2] = 1/(2*pi*sigma*sigma) * exp(-(i*i+j*j)/(2*sigma*sigma));
    }
    for (i=2; i<h-2; i++) for (j=2; j<w-2; j++) {
        res = 0.0;
        for (p=-2; p<3; p++) for (q=-2; q<3; q++) {
            res += G[p+2][q+2] * Pict[i+p][j+q].R;
        }
        Pict_new[i][j].R = res;
        Pict_new[i][j].G = res;
        Pict_new[i][j].B = res;
        Pict_new[i][j].alpha = 255;
    }
    return;
}

void Sobel_operator(int h, int w, pixel** Pict, pixel** Pict_border, int border_value)
{
    int i, j, Gx, Gy, G;
    for (i=1; i<h-1; i++) for (j=1; j<w-1; j++) {
        Gx = (Pict[i+1][j-1].R + 2*Pict[i+1][j].R + Pict[i+1][j+1].R) - (Pict[i-1][j-1].R + 2*Pict[i-1][j].R + Pict[i-1][j+1].R);
        Gy = (Pict[i-1][j+1].R + 2*Pict[i][j+1].R + Pict[i+1][j+1].R) - (Pict[i-1][j-1].R + 2*Pict[i][j-1].R + Pict[i+1][j-1].R);
        G = sqrt(Gx*Gx + Gy*Gy);
        if (G < border_value) G = 0;
        Pict_border[i][j].R = G;
        Pict_border[i][j].G = G;
        Pict_border[i][j].B = G;
        Pict_border[i][j].alpha = 255;
    }
    return;
}

void Canny_operator(int h, int w, pixel** Pict, pixel** Pict_border, int lower_border_value, int upper_border_value)
{
    int i, j, Gx, Gy, G;
    double theta;
    for (i=1; i<h-1; i++) for (j=1; j<w-1; j++) {
        Gx = (Pict[i+1][j-1].R + 2*Pict[i+1][j].R + Pict[i+1][j+1].R) - (Pict[i-1][j-1].R + 2*Pict[i-1][j].R + Pict[i-1][j+1].R);
        Gy = (Pict[i-1][j+1].R + 2*Pict[i][j+1].R + Pict[i+1][j+1].R) - (Pict[i-1][j-1].R + 2*Pict[i][j-1].R + Pict[i+1][j-1].R);
        G = sqrt(Gx*Gx + Gy*Gy);
        theta = atan2(Gy, Gx);
        if (theta < 0) theta += pi;
        if (theta < pi/8 || theta > 7*pi/8) if (G < Pict[i][j-1].R || G < Pict[i][j+1].R) G = 0;
        else if (theta >= pi/8 && theta < 3*pi/8) if (G < Pict[i-1][j+1].R || G < Pict[i+1][j-1].R) G = 0;
        else if (theta >= 3*pi/8 && theta < 5*pi/8) if (G < Pict[i-1][j].R || G < Pict[i+1][j].R) G = 0;
        else if (G < Pict[i+1][j+1].R || G < Pict[i-1][j-1].R) G = 0;
        if (G < lower_border_value) G = 0;
        Pict_border[i][j].R = G;
        Pict_border[i][j].G = G;
        Pict_border[i][j].B = G;
        Pict_border[i][j].alpha = 255;
    }
    return;
}

DUS* Make_Set(pixel* P)
{
    DUS* head = malloc(sizeof(DUS));
    Elem* el = malloc(sizeof(Elem));
    el->data = P;
    el->head = head;
    el->next = NULL;
    head->beg = el;
    head->end = el;
    head->color_status = 0;
    head->size = 1;
    return head;
}

DUS* Find_Set(Elem* el)
{
    return el->head;
}

DUS* Union_Set(Elem* el_1, Elem* el_2)
{
    DUS* s_1 = Find_Set(el_1), *s_2 = Find_Set(el_2);
    if (s_1 != s_2) {
        if (s_1->size < s_2->size) {
            DUS* s_3 = s_2;
            s_2 = s_1;
            s_1 = s_3;
        }
        s_1->end->next = s_2->beg;
        s_1->end = s_2->end;
        Elem* tmp;
        for (tmp=s_2->beg; tmp; tmp=tmp->next) {
            tmp->head = s_1;
        }
        s_1->size += s_2->size;
        free(s_2);
    }
    return s_1;
}

void Clear_Set(Elem** Verts, int h, int w)
{
    int i;
    Elem* tmp;
    for (i=0; i<h*w; i++) {
        if (Verts[i]->head) {
            tmp = Verts[i]->head->beg;
            free(tmp->head);
            while (tmp) {
                tmp->head = NULL;
                tmp = tmp->next;
            }
        }
        free(Verts[i]);
    }
    return;
}

void Make_DUS(int h, int w, pixel** Pict, int min_diff, int min_size, int max_size)
{
    int i, j, R, G, B;
    Elem** Verts = malloc(h*w*sizeof(Elem*)), *tmp = NULL;
    for (i=0; i<h; i++) for (j=0; j<w; j++) {
        Verts[i*w+j] = Make_Set(&Pict[i][j])->beg;
    }
    for (i=1; i<h; i++) for (j=1; j<w; j++) {
        if (abs(Pict[i][j].R - Pict[i][j-1].R) <= min_diff && abs(Pict[i][j].R - Pict[i-1][j].R) <= min_diff) {
            Union_Set(Verts[i*w+j], Verts[i*w+j-1]);
            Union_Set(Verts[i*w+j], Verts[i*w+j-w]);
        }
    }
    srand(time(NULL));
    int cur_color_status = 1;
    for (i=0; i<h*w; i++) if (Verts[i]->head->color_status == 0 && Verts[i]->head->size >= min_size && Verts[i]->head->size <= max_size) {
        R = ((double) rand()/RAND_MAX * 128.0) + 127;
        G = ((double) rand()/RAND_MAX * 128.0) + 127;
        B = ((double) rand()/RAND_MAX * 128.0) + 127;
        tmp = Verts[i]->head->beg;
        while (tmp) {
            tmp->data->R = R;
            tmp->data->G = G;
            tmp->data->B = B;
            tmp->data->alpha = 255;
            tmp = tmp->next;
        }
        Verts[i]->head->color_status = cur_color_status;
        cur_color_status++;
    }
    Clear_Set(Verts, h, w);
    free(Verts);
    return;
}

int main()
{
    int h=0, w=0;
    int border_value = 32; // граничное значение ЧБ цвета для отделения изображения от фона
    int min_diff = 10; // модуль максимальной разности цветов для соседних пикселей в одной компоненте связности
    int min_size = 50; // максимальный размер отбрасываемых при раскрашивании компонент связности (в пикселях)
    int max_size = 1000000; // максимальный размер раскрашиваемых компонент связности (в пикселях)
    double sigma = 0.84; // дисперсия размытия изображения под действием оператора Гаусса 5*5

	char *ins_file = "skull3.png"; // название файла с исходным изображением
	char *out_file = "skill3.png"; // название файла, в который нужно положить результат

    pixel** Pict = get_image(ins_file, &h, &w);

    image_to_BW(h, w, Pict, border_value);

    pixel** Pict_new = create_image(h, w, NULL);
    Gauss_operator(h, w, Pict, Pict_new, sigma);

    pixel** Pict_border = create_image(h, w, NULL);
//    Canny_operator(h, w, Pict_new, Pict_border, border_value, 100);
    Sobel_operator(h,w, Pict_new, Pict_border, border_value);

    Make_DUS(h, w, Pict_border, min_diff, min_size, max_size);

    put_image(h, w, Pict_border, out_file);
    printf("The answer image is ready! It is inside this folder nicknamed '%s'.", out_file);

    clear_image(h, w, Pict);
    clear_image(h, w, Pict_new);
    clear_image(h, w, Pict_border);
	return 0;
}
