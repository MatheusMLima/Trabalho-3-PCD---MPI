#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/time.h>

#define N 2048
#define NUM_GER 2000

int NUM_THREADS;

int vivos = 0; // global pra não ter que ficar passando por parametro esse treco

// Retorna a diferença de tempo entre "start" e "end"
float get_tempo(struct timespec *start, struct timespec *end)
{
    float time_taken;
    time_taken = (end->tv_sec - start->tv_sec) + 1e-9 * (end->tv_nsec - start->tv_nsec);
    return time_taken;
}

void iniciaMatriz(float ***tabuleiro)
{
    *tabuleiro = (float **)calloc(N, sizeof(float *));
    int i;
    for (i = 0; i < N; i++)
    {
        (*tabuleiro)[i] = (float *)calloc(N, sizeof(float));
    }
}

void iniciaGlider(float ***tabuleiro)
{
    int lin = 1, col = 1;
    (*tabuleiro)[lin][col + 1] = 1.0;
    (*tabuleiro)[lin + 1][col + 2] = 1.0;
    (*tabuleiro)[lin + 2][col] = 1.0;
    (*tabuleiro)[lin + 2][col + 1] = 1.0;
    (*tabuleiro)[lin + 2][col + 2] = 1.0;
}

void iniciaPentomino(float ***tabuleiro)
{
    int lin = 10, col = 30;
    (*tabuleiro)[lin][col + 1] = 1.0;
    (*tabuleiro)[lin][col + 2] = 1.0;
    (*tabuleiro)[lin + 1][col] = 1.0;
    (*tabuleiro)[lin + 1][col + 1] = 1.0;
    (*tabuleiro)[lin + 2][col + 1] = 1.0;
}

void imprimeMatriz(float **tabuleiro)
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            printf("%.1f | ", tabuleiro[i][j]);
        }
        printf("\n");
    }
}

int getNeighbors(float **tabuleiro, int lin, int col)
{

    int cima, baixo, esquerda, direita;

    cima = lin - 1;
    baixo = lin + 1;
    esquerda = col - 1;
    direita = col + 1;

    if (cima < 0)
        cima = N - 1;

    if (baixo >= N)
        baixo = 0;

    if (esquerda < 0)
        esquerda = N - 1;

    if (direita >= N)
        direita = 0;

    return tabuleiro[cima][esquerda] + tabuleiro[cima][col] + tabuleiro[cima][direita] + tabuleiro[lin][esquerda] + tabuleiro[lin][direita] + tabuleiro[baixo][esquerda] + tabuleiro[baixo][col] + tabuleiro[baixo][direita];
}

int AtualizaGrid(float **grid, float **newgrid, int rank)
{
    int i, j, vizinhos_vivos, total_vivos = 0;
    for (i = (N / NUM_THREADS) * rank; i < (N / NUM_THREADS) * (rank + 1); i++)
    {
        for (j = 0; j < N; j++)
        {
            vizinhos_vivos = getNeighbors(grid, i, j);
            if (grid[i][j] == 1 && (vizinhos_vivos == 2 || vizinhos_vivos == 3))
                newgrid[i][j] = 1;
            else if (grid[i][j] == 0 && vizinhos_vivos == 3)
                newgrid[i][j] = 1;
            else
                newgrid[i][j] = grid[i][j] * 0;

            total_vivos += newgrid[i][j];
        }
    }
    return total_vivos;
}

int envio_msg(int i, int j)
{
    int modo = i % j;
    while (modo < 0)
    {
        modo = (j + modo) % j;
    }
    return modo;
}

void mensagem(float **grid, int rank, int modo, int geracao)
{
    MPI_Request request_s, request_r;
    MPI_Status status_s, status_r;
    int position_r, position_s, rank_s, rank_r, neighbor_process;
    if (modo == 0)
    {
        position_s = (N / NUM_THREADS) * (rank + 1) - 1;
        position_r = envio_msg((N / NUM_THREADS) * rank - 1, N);
        rank_s = envio_msg(rank + 1, NUM_THREADS);
        rank_r = envio_msg(rank - 1, NUM_THREADS);
    }
    else if (modo == 1)
    {
        position_s = (N / NUM_THREADS) * rank;
        position_r = envio_msg((N / NUM_THREADS) * (rank + 1), N);
        rank_s = envio_msg(rank - 1, NUM_THREADS);
        rank_r = envio_msg(rank + 1, NUM_THREADS);
    }

    neighbor_process = geracao * 2 + 1 + modo;

    MPI_Isend(grid[position_s], N, MPI_INT, rank_s, neighbor_process, MPI_COMM_WORLD, &request_s);
    MPI_Irecv(grid[position_r], N, MPI_INT, rank_r, neighbor_process, MPI_COMM_WORLD, &request_r);
    MPI_Wait(&request_s, &status_s);
    MPI_Wait(&request_r, &status_r);
}

void start_game(float **grid, float **newgrid, int rank, int local_alives, int global_alives, double start)
{
    double end;
    if (rank == 0)
    {
        start = MPI_Wtime();
    }

    for (int i = 0; i < NUM_GER; i++)
    {
        local_alives = AtualizaGrid(grid, newgrid, rank);
        MPI_Reduce(&local_alives, &global_alives, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank == 0)
        {
            printf("Geracao %d: %d\n", i + 1, global_alives);
        }

        float **j = grid;
        grid = newgrid;
        newgrid = j;
        if (NUM_THREADS > 1)
        {
            mensagem(grid, rank, 0, i);
            mensagem(grid, rank, 1, i);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
    {
        printf("Células vivas na última geração: %d\n", global_alives);
        end = MPI_Wtime();
        printf("Tempo decorrido na sessão paralela: %0.3f\n", end - start);
    }

    for (int i = 0; i < N; i++)
    {
        free(grid[i]);
        free(newgrid[i]);
    }
}

int main()
{
    MPI_Init(NULL, NULL);
    int rank, local_alives, global_alives;
    double start = 0.0;

    float **grid, **newGrid;

    // start = MPI_Wtime();

    iniciaMatriz(&grid);
    iniciaMatriz(&newGrid);

    MPI_Comm_size(MPI_COMM_WORLD, &NUM_THREADS);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    iniciaGlider(&grid);
    iniciaPentomino(&grid);

    start_game(grid, newGrid, rank, local_alives, global_alives, start);

    free(grid);
    free(newGrid);

    MPI_Finalize();
    return 0;
}