#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <mpi.h>

#define MASTER 0

typedef struct {
    int alive;
    int row;
    int col;
} Invader;

typedef struct {
    int active;
    int col;
    int fire_tick;
    int from_row;
} Cannonball;

void print_game(int n, int m, Invader* invaders, int player_col, int tick);

int main(int argc, char *argv[]) {
    int rank, size;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (argc != 3) {
        if (rank == MASTER) {
            printf("Usage: %s <rows> <cols>\n", argv[0]);
            printf("Example: mpirun -np 7 --oversubscribe %s 3 2\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }
    
    int n = atoi(argv[1]);  // rows
    int m = atoi(argv[2]);  // cols
    
    if (size != n * m + 1) {
        if (rank == MASTER) {
            printf("ERROR: Need exactly %d processes (1 player + %d invaders)\n", 
                   n * m + 1, n * m);
            printf("You provided %d processes\n", size);
        }
        MPI_Finalize();
        return 1;
    }
    
    // MASTER PROCESS (Player)
    if (rank == MASTER) {
        printf("========================================\n");
        printf("   SPACE INVADERS MPI SIMULATION\n");
        printf("========================================\n");
        printf("Grid: %d rows x %d columns\n", n, m);
        printf("Total invaders: %d\n", n * m);
        printf("Player starts at column 0\n");
        printf("\nTravel time:\n");
        printf("  Player -> Bottom row: 2 seconds\n");
        printf("  Player -> Top row: %d seconds\n", n + 1);
        printf("  Invaders shoot every 4 seconds (10%% chance)\n\n");
        
        int player_col = 0;
        int tick = 0;
        int game_over = 0;
        int player_alive = 1;
        
        Invader invaders[n * m];
        for (int i = 0; i < n * m; i++) {
            invaders[i].alive = 1;
            invaders[i].row = i / m;
            invaders[i].col = i % m;
        }
        
        Cannonball player_cannonballs[100];
        int num_player_cannonballs = 0;
        
        Cannonball enemy_cannonballs[100];
        int num_enemy_cannonballs = 0;
        
        srand(time(NULL));
        
        while (!game_over) {
            tick++;
            printf("\n========== TICK %d ==========\n", tick);
            
            // 1. PLAYER MOVEMENT
            int move = rand() % 3 - 1;
            int new_col = player_col + move;
            
            if (new_col < 0 || new_col >= m) {
                new_col = player_col;
            }
            player_col = new_col;
            
            printf("Player at column %d\n", player_col);
            
            // 2. PLAYER FIRES
            Cannonball new_ball;
            new_ball.active = 1;
            new_ball.col = player_col;
            new_ball.fire_tick = tick;
            new_ball.from_row = -1;
            player_cannonballs[num_player_cannonballs++] = new_ball;
            
            // 3. RECEIVE STATUS FROM ALL INVADERS
            for (int i = 1; i <= n * m; i++) {
                int status;
                MPI_Recv(&status, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                invaders[i - 1].alive = status;
            }
            
            // 4. RECEIVE FIRE INFO
            for (int i = 1; i <= n * m; i++) {
                int fired_tick;
                MPI_Recv(&fired_tick, 1, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                if (fired_tick > 0) {
                    int inv_row = invaders[i - 1].row;
                    int inv_col = invaders[i - 1].col;
                    printf("  Invader at (%d,%d) fired at tick %d\n", inv_row, inv_col, fired_tick);
                    
                    Cannonball eb;
                    eb.active = 1;
                    eb.col = inv_col;
                    eb.fire_tick = fired_tick;
                    eb.from_row = inv_row;
                    enemy_cannonballs[num_enemy_cannonballs++] = eb;
                }
            }
            
            // 5. CHECK PLAYER CANNONBALLS FOR HITS
            for (int i = 0; i < num_player_cannonballs; i++) {
                if (!player_cannonballs[i].active) continue;
                
                int time_traveled = tick - player_cannonballs[i].fire_tick;
                
                for (int row = n - 1; row >= 0; row--) {
                    int time_needed = 2 + (n - 1 - row);
                    
                    if (time_traveled == time_needed) {
                        int target_index = row * m + player_cannonballs[i].col;
                        
                        if (invaders[target_index].alive) {
                            int is_bottom_most = 1;
                            for (int r = row + 1; r < n; r++) {
                                int check_index = r * m + player_cannonballs[i].col;
                                if (invaders[check_index].alive) {
                                    is_bottom_most = 0;
                                    break;
                                }
                            }
                            
                            if (is_bottom_most) {
                                printf("  >>> HIT! Player destroyed invader at (%d,%d) [travel: %d sec]\n", 
                                       row, player_cannonballs[i].col, time_needed);
                                
                                int target_rank = target_index + 1;
                                int kill = 1;
                                MPI_Send(&kill, 1, MPI_INT, target_rank, 2, MPI_COMM_WORLD);
                                
                                invaders[target_index].alive = 0;
                                player_cannonballs[i].active = 0;
                                break;
                            }
                        }
                    }
                }
            }
            
            // 6. CHECK ENEMY CANNONBALLS
            for (int i = 0; i < num_enemy_cannonballs; i++) {
                if (!enemy_cannonballs[i].active) continue;
                
                int time_traveled = tick - enemy_cannonballs[i].fire_tick;
                int inv_row = enemy_cannonballs[i].from_row;
                int inv_col = enemy_cannonballs[i].col;
                int inv_index = inv_row * m + inv_col;
                
                if (!invaders[inv_index].alive) continue;
                
                int is_bottom_most = 1;
                for (int r = inv_row + 1; r < n; r++) {
                    int check_index = r * m + inv_col;
                    if (invaders[check_index].alive) {
                        is_bottom_most = 0;
                        break;
                    }
                }
                
                if (!is_bottom_most) continue;
                
                int time_needed = 2 + (n - 1 - inv_row);
                
                if (time_traveled == time_needed && player_col == inv_col) {
                    printf("  !!! PLAYER HIT by invader at (%d,%d) [travel: %d sec] !!!\n", 
                           inv_row, inv_col, time_needed);
                    player_alive = 0;
                    game_over = 1;
                    enemy_cannonballs[i].active = 0;
                    break;
                }
            }
            
            // 7. COUNT ALIVE
            int alive_count = 0;
            for (int i = 0; i < n * m; i++) {
                if (invaders[i].alive) alive_count++;
            }
            printf("Invaders remaining: %d/%d\n", alive_count, n * m);
            
            if (alive_count == 0) {
                printf("\n*** VICTORY! All invaders destroyed! ***\n");
                game_over = 1;
            }
            
            // 8. PRINT GAME
            print_game(n, m, invaders, player_col, tick);
            
            // 9. SYNCHRONIZE
            MPI_Barrier(MPI_COMM_WORLD);
            
            if (tick >= 100) {
                game_over = 1;
            }
            
            sleep(1);
        }
        
        // Stop all invaders
        for (int i = 1; i <= n * m; i++) {
            int stop = 1;
            MPI_Send(&stop, 1, MPI_INT, i, 99, MPI_COMM_WORLD);
        }
        
        printf("\n========================================\n");
        if (!player_alive) {
            printf("         GAME OVER - DEFEAT\n");
        } else {
            printf("        GAME OVER - VICTORY\n");
        }
        printf("========================================\n");
        
    } else {
        // INVADER PROCESS - SIMPLIFIED WITHOUT MPI_Allreduce
        int my_row = (rank - 1) / m;
        int my_col = (rank - 1) % m;
        int alive = 1;
        int my_tick = 0;
        
        srand(time(NULL) + rank);
        
        while (1) {
            my_tick++;
            
            // 1. Send status
            MPI_Send(&alive, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
            
            // 2. Decide if firing
            // Note: We let master check if bottom-most, so any alive invader can attempt to fire
            int fired_tick = 0;
            if (alive && my_tick % 4 == 0) {
                if ((rand() % 100) < 10) {
                    fired_tick = my_tick;
                }
            }
            MPI_Send(&fired_tick, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD);
            
            // 3. Check if killed
            int flag;
            MPI_Status status;
            MPI_Iprobe(MASTER, 2, MPI_COMM_WORLD, &flag, &status);
            if (flag) {
                int kill;
                MPI_Recv(&kill, 1, MPI_INT, MASTER, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                alive = 0;
            }
            
            // 4. Check for stop
            MPI_Iprobe(MASTER, 99, MPI_COMM_WORLD, &flag, &status);
            if (flag) {
                int stop;
                MPI_Recv(&stop, 1, MPI_INT, MASTER, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                break;
            }
            
            // 5. Synchronize
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    
    MPI_Finalize();
    return 0;
}

void print_game(int n, int m, Invader* invaders, int player_col, int tick) {
    printf("\n--- Game Board ---\n");
    
    for (int row = 0; row < n; row++) {
        printf("Row %d: ", row);
        for (int col = 0; col < m; col++) {
            int index = row * m + col;
            if (invaders[index].alive) {
                printf("[👾 ] ");
            } else {
                printf("[ ] ");
            }
        }
        printf("\n");
    }
    
    printf("\n       ");
    for (int col = 0; col < m; col++) {
        if (col == player_col) {
            printf(" ^  ");
        } else {
            printf("    ");
        }
    }
    printf("\n");
    printf("Player:");
    for (int col = 0; col < m; col++) {
        if (col == player_col) {
            printf("[P] ");
        } else {
            printf("    ");
        }
    }
    printf("\n------------------\n");
}
