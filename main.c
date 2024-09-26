#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define TAM 48

int ROWS;
int COLS;

// u(t) : N = β0 + β1t
typedef struct regressaoLinear {
  double B0;
  double B1;
} RegressaoLinear;

// u(t) : N = β0 + β1t + β2t
typedef struct regressaoQuadratica
{
  double B0;
  double B1;
  double B2;
} RegressaoQuadratica;

// u(t) : N = β0e^(β1t)
typedef struct regressaoNaoLinear
{
  double B0;
  double B1;
} RegressaoNaoLinear;

typedef struct dados {
  RegressaoLinear* rLinear;
  RegressaoQuadratica* rQuadratica;
  RegressaoNaoLinear* rNaoLinear;
} Dados;
void print_matrix(double matrix[ROWS][COLS]);
void gaussElimination(double M[ROWS][COLS], Dados* dados, int tipo);
double regressaoLinear(double matrixDeSomatorio[ROWS][COLS], double vetorX[TAM], double vetorY[TAM], double somatorioDeYAoQuadrado, double somatorioDeY);
double regressaoQuadratica(double matrixDeSomatorio[ROWS][COLS], double vetorX[TAM], double vetorY[TAM], double vetorXAoQuadrado[TAM], double somatorioDeYAoQuadrado, double somatorioDeY);
double regressaoNaoLinear(double matrixDeSomatorio[ROWS][COLS], double vetorX[TAM], double vetorY[TAM], double somatorioDeYAoQuadrado, double somatorioDeY);

int main(char argc, char argv[]) {
  FILE* arquivo = fopen("dados.txt", "r");

  double* vetorX = malloc(sizeof(double) * TAM);
  double* xAoQuadrado = malloc(sizeof(double) * TAM);
  double* xAoCubo = malloc(sizeof(double) * TAM);
  double* xAQuarta = malloc(sizeof(double) * TAM);

  double* vetorY = malloc(sizeof(double) * TAM);
  double* yAoQuadrado = malloc(sizeof(double) * TAM);
  double* xVezesY = malloc(sizeof(double) * TAM);
  double* LogNaturalDeY = malloc(sizeof(double) * TAM);
  double* XVezesLogNaturalDeY = malloc(sizeof(double) * TAM);

  double somatorioDeX = 0.0;
  double somatorioDeXAoQuadrado = 0.0;
  double somatorioDeXAoCubo = 0.0;
  double somatorioDeXAQuarta = 0.0;

  double somatorioDeY = 0.0;
  double somatorioDeXVezesY = 0.0;
  double somatorioDeYVezesXAoQuadrado = 0.0;
  double somatorioDeYAoQuadrado = 0.0;
  double somatorioDeLogNaturalDeY = 0.0;
  double somatorioDeXVezesLogNaturalDeY = 0.0;

  if (arquivo == NULL){
    printf("Erro de alocacao ou Falha no arquivo \n");
    exit(1);
  }

  while (feof(arquivo) == 0){
    for (int i = 0; i < TAM; i++){
      fscanf(arquivo, "%lf %lf\n", &vetorX[i], &vetorY[i]);
      xAoQuadrado[i] = vetorX[i] * vetorX[i];
      xAoCubo[i] = pow(vetorX[i], 3);
      xAQuarta[i] = pow(vetorX[i], 4);
      xVezesY[i] = vetorX[i] * vetorY[i];
      yAoQuadrado[i] = pow(vetorY[i], 2);
      LogNaturalDeY[i] = log(vetorY[i]);
      XVezesLogNaturalDeY[i] = vetorX[i] * log(vetorY[i]);
    }
  }

  if (arquivo == NULL){
    printf("Erro de alocacao ou Falha no arquivo \n");
    exit(1);
  }

  for (int i = 0; i < TAM; i++){
    somatorioDeX += vetorX[i];
    somatorioDeXAoQuadrado += xAoQuadrado[i];
    somatorioDeXAoCubo += xAoCubo[i];
    somatorioDeXAQuarta += xAQuarta[i];

    somatorioDeY += vetorY[i];
    somatorioDeYAoQuadrado += yAoQuadrado[i];
    somatorioDeXVezesY += xVezesY[i];
    somatorioDeYVezesXAoQuadrado += vetorY[i] * xAoQuadrado[i];
    somatorioDeLogNaturalDeY += LogNaturalDeY[i];
    somatorioDeXVezesLogNaturalDeY += XVezesLogNaturalDeY[i];
  }

  ROWS = 2;
  COLS = 3;

  double matrizRegressaoLinear[2][3] = {
      {TAM, somatorioDeX, somatorioDeY},
      {somatorioDeX, somatorioDeXAoQuadrado, somatorioDeXVezesY},
  };

  double cDeterminacaoRegressaoLinear = regressaoLinear(matrizRegressaoLinear, vetorX, vetorY, somatorioDeYAoQuadrado, somatorioDeY);

  ROWS = 3;
  COLS = 4;

  double matrizRegressaoQuadratica[3][4] = {
      {TAM, somatorioDeX, somatorioDeXAoQuadrado, somatorioDeY},
      {somatorioDeX, somatorioDeXAoQuadrado, somatorioDeXAoCubo, somatorioDeXVezesY},
      {somatorioDeXAoQuadrado, somatorioDeXAoCubo, somatorioDeXAQuarta, somatorioDeYVezesXAoQuadrado},
  };

  double cDeterminacaoRegressaoQuadratica = regressaoQuadratica(matrizRegressaoQuadratica, vetorX, vetorY, xAoQuadrado, somatorioDeYAoQuadrado, somatorioDeY);

  ROWS = 2;
  COLS = 3;

  double matrizRegressaoNaoLinear[2][3] = {
      {TAM, somatorioDeX, somatorioDeLogNaturalDeY},
      {somatorioDeX, somatorioDeXAoQuadrado, somatorioDeXVezesLogNaturalDeY},
  };

  double cDeterminacaoRegressaoNaoLinear = regressaoNaoLinear(matrizRegressaoNaoLinear, vetorX, vetorY, somatorioDeYAoQuadrado, somatorioDeY);

  printf("\n");

  if(cDeterminacaoRegressaoLinear > cDeterminacaoRegressaoQuadratica && cDeterminacaoRegressaoLinear > cDeterminacaoRegressaoNaoLinear) {
    printf("Regressao Linear melhor se ajustou aos dados\nr^2 = %lf", cDeterminacaoRegressaoLinear);
  }
  else if (cDeterminacaoRegressaoQuadratica > cDeterminacaoRegressaoNaoLinear){
    printf("Regressao Quadratica melhor se ajustou aos dados\nr^2 = %lf", cDeterminacaoRegressaoQuadratica);
  }
  else{
    printf("Regressao Não Linear melhor se ajustou aos dados\nr^2 = %lf", cDeterminacaoRegressaoNaoLinear);
  }

  fclose(arquivo);
  return 0;
}
void print_matrix(double matrix[ROWS][COLS]){
  for (int i = 0; i < ROWS; i++){
    for (int j = 0; j < COLS; j++){
      printf("%lf ", matrix[i][j]);
    }
    printf("\n");
  }
}
//tipo == 0: Regressão Linear | tipo == 1: Regressão Quadrática | tipo == 2: Regressão Não Linear
void gaussElimination(double M[ROWS][COLS], Dados* dados, int tipo){
  double M2[ROWS][COLS];

  memcpy(M2, M, ROWS * COLS * sizeof(double));

  for (int j = 0; j < COLS - 2; j++){
    for (int i = j; i < ROWS; i++){
      if (M2[i][j] != 0){
        if (i != j){
          for (int k = 0; k < COLS; k++){
            double temp = M2[i][k];
            M2[i][k] = M2[j][k];
            M2[j][k] = temp;
          }
        }
        for (int m = j + 1; m < ROWS; m++){
          double a = -M2[m][j] / M2[j][j];
          for (int n = j; n < COLS; n++){
            M2[m][n] += a * M2[j][n];
          }
        }
        printf("\n");

        break;
      }
    }
  }

  double x[ROWS];

  x[ROWS - 1] = M2[ROWS - 1][ROWS] / M2[ROWS - 1][ROWS - 1];

  for (int i = ROWS - 2; i >= 0; i--){
    double sum = 0.0;
    for (int j = i + 1; j < ROWS; j++){
      sum = sum + M2[i][j] * x[j];
    }
    x[i] = (M2[i][ROWS] - sum) / M2[i][i];
  }

  if(tipo == 0) {
    dados->rLinear->B0 = x[0];
    dados->rLinear->B1 = x[1];
  }
  else if(tipo == 1) {
    dados->rQuadratica->B0 = x[0];
    dados->rQuadratica->B1 = x[1];
    dados->rQuadratica->B2 = x[2];
  }
  else if(tipo == 2) {
    dados->rNaoLinear->B0 = exp(x[0]);     //Esse é o único caso onde Bi não é, diretamente, x[i]. B0 na verdade é log de x[0] na base e (número de Euler), ou seja B0 = exp(x[0]), portanto precisamos fazer uma substituição
    dados->rNaoLinear->B1 = x[1];
  }
}

double regressaoLinear(double matrixDeSomatorio[ROWS][COLS], double vetorX[TAM], double vetorY[TAM], double somatorioDeYAoQuadrado, double somatorioDeY){
  Dados dados;
  dados.rLinear = malloc(sizeof(RegressaoLinear));

  double *vetorUi = malloc(sizeof(double) * TAM);
  double *vetorDi = malloc(sizeof(double) * TAM);

  double somatorioDiAoQuadrado = 0.0;
  double rAoQuadrado = 0.0;

  gaussElimination(matrixDeSomatorio, &dados, 0);
  printf("Regressao Linear\n");
  printf("b0: %lf\nb1: %lf\n", dados.rLinear->B0, dados.rLinear->B1);

  for (int i = 0; i < TAM; i++){
    vetorUi[i] = dados.rLinear->B0 + dados.rLinear->B1 * vetorX[i];
    vetorDi[i] = vetorY[i] - vetorUi[i];
    somatorioDiAoQuadrado += pow(vetorDi[i], 2);
  }

  double tamanho = TAM;
  double temp = (somatorioDeYAoQuadrado - (pow(somatorioDeY, 2)) / tamanho);

  rAoQuadrado = 1 - (somatorioDiAoQuadrado / temp);
  printf("Coeficiente de Determinacao: %lf\n", rAoQuadrado);
  return rAoQuadrado;
}

double regressaoQuadratica(double matrixDeSomatorio[ROWS][COLS], double vetorX[TAM], double vetorY[TAM], double vetorXAoQuadrado[TAM], double somatorioDeYAoQuadrado, double somatorioDeY){
  Dados dados;
  dados.rQuadratica = malloc(sizeof(RegressaoQuadratica));

  double *vetorUi = malloc(sizeof(double) * TAM);
  double *vetorDi = malloc(sizeof(double) * TAM);
  double somatorioVetorDiAoQuadrado = 0.0;

  double rAoQuadrado = 0.0;

  gaussElimination(matrixDeSomatorio, &dados, 1);
  printf("Regressao Quadratica\n");
  printf("b0: %lf\nb1: %lf\nb2: %lf\n", dados.rQuadratica->B0, dados.rQuadratica->B1, dados.rQuadratica->B2);

  for (int i = 0; i < TAM; i++){
    vetorUi[i] = dados.rQuadratica->B0 + (dados.rQuadratica->B1 * vetorX[i]) + (dados.rQuadratica->B2 * vetorXAoQuadrado[i]);
    vetorDi[i] = vetorY[i] - vetorUi[i];
    somatorioVetorDiAoQuadrado += pow(vetorDi[i], 2);
  }

  double tamanho = TAM;
  double temp = (somatorioDeYAoQuadrado - (pow(somatorioDeY, 2)) / tamanho);

  rAoQuadrado = 1 - (somatorioVetorDiAoQuadrado / temp);
  printf("Coeficiente de Determinacao: %lf\n", rAoQuadrado);
  return rAoQuadrado;
}

double regressaoNaoLinear(double matrixDeSomatorio[ROWS][COLS], double vetorX[TAM], double vetorY[TAM], double somatorioDeYAoQuadrado, double somatorioDeY){
  Dados dados;
  dados.rNaoLinear = malloc(sizeof(RegressaoNaoLinear));

  double *vetorUi = malloc(sizeof(double) * TAM);
  double *vetorDi = malloc(sizeof(double) * TAM);

  double somatorioDiAoQuadrado = 0.0;
  double rAoQuadrado = 0.0;

  gaussElimination(matrixDeSomatorio, &dados, 2);
  printf("Regressao Nao Linear\n");
  printf("b0: %lf\nb1: %lf\n", dados.rNaoLinear->B0, dados.rNaoLinear->B1);
  for (int i = 0; i < TAM; i++){
    vetorUi[i] = dados.rNaoLinear->B0 * exp(dados.rNaoLinear->B1 * vetorX[i]);
    vetorDi[i] = vetorY[i] - vetorUi[i];
    somatorioDiAoQuadrado += pow(vetorDi[i], 2);
  }

  double tamanho = TAM;
  double temp = (somatorioDeYAoQuadrado - (pow(somatorioDeY, 2)) / tamanho);

  rAoQuadrado = 1 - (somatorioDiAoQuadrado / temp);
  printf("Coeficiente de Determinacao: %lf\n", rAoQuadrado);
  return rAoQuadrado;
}
