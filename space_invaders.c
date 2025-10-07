#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <mpi.h>

#define NUMPROCS_ARG 0
#define ROW_ARG 1
#define COL_ARG 2


#define MASTER 0
#define STATUS_TAG 0
#define INVADER_SHOT_TAG 1
#define PLAYER_SHOT_TAG 2
#define KILL_TAG 3


//structure to track invader state
typedef struct {
    int alive;
    int row;
    int col;
} Invader;

//structure for cannonballs
typedef struct {
    int active;
    int col;
    int fire_tick;
    int from_row;  // -1 means from player
} Cannonball;

//structure for players
typedef struct {
    int alive;
    int col;
} PlayerStruct;

void print_game(int n, int m, Invader* invaders, int player_col, Cannonball* player_balls, 
                int num_pballs, Cannonball* enemy_balls, int num_eballs);

int main(int argc, char *argv[]) {
    int rank, size;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    if (argc != 3) {
        if (rank == MASTER) {
            printf("Usage: %s <rows> <cols>\n", argv[NUMPROCS_ARG]);
            printf("Example: mpirun -np 7 --oversubscribe %s 3 2\n", argv[NUMPROCS_ARG]);
        }
        MPI_Finalize();
        return 1;
    }
    
    int n = atoi(argv[ROW_ARG]);  // rows
    int m = atoi(argv[COL_ARG]);  // cols
    

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
        printf("Player starts at column 0\n\n");
        
        int player_col = 0;
        int tick = 0;
        int game_over = 0;
        int player_alive = 1;
        
        // Track  invaders
        Invader invaders[n * m];
        for (int i = 0; i < n * m; i++) {
            invaders[i].alive = 1;
            invaders[i].row = i / m;
            invaders[i].col = i % m;
        }
        
        // Track cannonballs
        Cannonball player_cannonballs[100];
        int num_player_cannonballs = 0;
        
        Cannonball enemy_cannonballs[100];
        int num_enemy_cannonballs = 0;
        
        srand(time(NULL));
        
        // Main game loop
        while (!game_over) {
            tick++;
            printf("\n========== TICK %d ==========\n", tick);
            
            // 1. PLAYER MOVEMENT (random for simulation)
            int move = rand() % 3 - 1;  // -1, 0, or 1
            player_col = player_col + move;
            
            // Keep player in bounds
            if (player_col < 0) player_col = 0;
            if (player_col >= m) player_col = m - 1;
            
            printf("Player moved to column %d\n", player_col);
            
            // 2. PLAYER FIRES CANNONBALL
            Cannonball new_ball;
            new_ball.active = 1;
            new_ball.col = player_col;
            new_ball.fire_tick = tick;
            new_ball.from_row = -1;  // -1 means from player
            player_cannonballs[num_player_cannonballs] = new_ball;
            num_player_cannonballs++;
            
            // 3. GET STATUS FROM ALL INVADERS
            for (int i = 1; i <= n * m; i++) {
                int status;
                MPI_Recv(&status, 1, MPI_INT, i, STATUS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                invaders[i - 1].alive = status;
            }
            
            // 4. CHECK IF INVADERS FIRED
            for (int i = 1; i <= n * m; i++) {
                int fired;
                MPI_Recv(&fired, 1, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                if (fired > 0) {
                    int inv_row = invaders[i - 1].row;
                    int inv_col = invaders[i - 1].col;
                    
                    printf("Invader at (%d,%d) fired!\n", inv_row, inv_col);
                    
                    Cannonball eb;
                    eb.active = 1;
                    eb.col = inv_col;
                    eb.fire_tick = fired;
                    eb.from_row = inv_row;
                    enemy_cannonballs[num_enemy_cannonballs] = eb;
                    num_enemy_cannonballs++;
                }
            }
            
            // 5. CHECK PLAYER CANNONBALLS FOR HITS
            for (int i = 0; i < num_player_cannonballs; i++) {
                if (!player_cannonballs[i].active) continue;
                
                int time_traveled = tick - player_cannonballs[i].fire_tick;
                
                // Check each row from bottom to top
                for (int row = n - 1; row >= 0; row--) {
                    // Calculate time needed to reach this row
                    int time_needed = 2 + (n - 1 - row);
                    
                    if (time_traveled == time_needed) {
                        int target_index = row * m + player_cannonballs[i].col;
                        
                        // Check if this invader is alive
                        if (invaders[target_index].alive) {
                            // Check if it's the bottom-most alive invader in this column
                            int is_bottom_most = 1;
                            for (int r = row + 1; r < n; r++) {
                                int check_index = r * m + player_cannonballs[i].col;
                                if (invaders[check_index].alive) {
                                    is_bottom_most = 0;
                                    break;
                                }
                            }
                            
                            if (is_bottom_most) {
                                printf(">>> HIT! Player destroyed invader at (%d,%d)\n", 
                                       row, player_cannonballs[i].col);
                                
                                // Tell invader it's dead
                                int rank = target_index + 1;
                                int kill = 1;
                                MPI_Send(&kill, 1, MPI_INT, rank, PLAYER_SHOT_TAG, MPI_COMM_WORLD);
                                
                                //invaders[target_index].alive = 0;
                                player_cannonballs[i].active = 0;
                                break;
                            }
                        }
                    }
                }
            }
            
            // 6. CHECK ENEMY CANNONBALLS HITTING PLAYER
            for (int i = 0; i < num_enemy_cannonballs; i++) {
                if (!enemy_cannonballs[i].active) continue;
                
                int time_traveled = tick - enemy_cannonballs[i].fire_tick;
                int inv_row = enemy_cannonballs[i].from_row;
                int inv_col = enemy_cannonballs[i].col;
                int inv_index = inv_row * m + inv_col;
                
                // Only count if invader is still alive
                if (!invaders[inv_index].alive) continue;
                
                // Check if this invader is bottom-most in its column
                int is_bottom_most = 1;
                for (int r = inv_row + 1; r < n; r++) {
                    int check_index = r * m + inv_col;
                    if (invaders[check_index].alive) {
                        is_bottom_most = 0;
                        break;
                    }
                }
                
                if (!is_bottom_most) continue;
                
                // Calculate time needed
                int time_needed = 2 + (n - 1 - inv_row);
                
                if (time_traveled == time_needed) {
                    // Check if player is in same column
                    if (player_col == inv_col) {
                        printf("\n!!! PLAYER HIT by invader at (%d,%d) !!!\n", inv_row, inv_col);
                        player_alive = 0;
                        game_over = 1;
                        enemy_cannonballs[i].active = 0;
                        break;
                    }
                }
            }
            
            // 7. CHECK WIN CONDITION
            int alive_count = 0;
            for (int i = 0; i < n * m; i++) {
                if (invaders[i].alive) alive_count++;
            }
            
            printf("Invaders remaining: %d\n", alive_count);
            
            if (alive_count == 0) {
                printf("\n*** VICTORY! All invaders destroyed! ***\n");
                game_over = 1;
            }
            

            print_game(n, m, invaders, player_col, player_cannonballs, 
                      num_player_cannonballs, enemy_cannonballs, num_enemy_cannonballs);
            
            // 9. SYNCHRONIZE ALL PROCESSES
            MPI_Barrier(MPI_COMM_WORLD);
            
            // Safety limit
            if (tick >= 100) {
                printf("\nReached tick limit\n");
                game_over = 1;
            }
            
            sleep(1);  // 1 second per tick
        }
        
        // Game over
        printf("\nGame Over! Sending stop signal to all invaders...\n");
        for (int i = 1; i <= n * m; i++) {
            int stop = 1;
            MPI_Send(&stop, 1, MPI_INT, i, 99, MPI_COMM_WORLD);
        }
        
        if (!player_alive) {
            printf("\n=== DEFEAT ===\n");
        } else {
            printf("\n=== VICTORY ===\n");
        }
        
    } else {
        // INVADER PROCESS
        int my_row = (rank - 1) / m;
        int my_col = (rank - 1) % m;
        int left = rank - 1;
        int right = rank + 1;
        int alive = 1;
        int my_tick = 0;
        
        srand(time(NULL) + rank);
        if (my_col == 0){
            left = MPI_PROC_NULL;
        }else if(my_col == m-1){
            right = MPI_PROC_NULL;
        }

        while (1) {
            my_tick++;
            
            // 1. Send status to master
            MPI_Send(&alive, 1, MPI_INT, MASTER, STATUS_TAG, MPI_COMM_WORLD);
            
            // 2. Decide if firing (10% chance every 4 ticks)
            int fired = 0;
            if (alive && my_tick % 4 == 0) {
                int random = rand() % 100;
                if (random < 10) {  // 10% probability
                    fired = my_tick;
                }
            }
            MPI_Send(&fired, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD);
            

            int flag;
            MPI_Status status;
            //Invader has been shot by the player
            MPI_Iprobe(MASTER, PLAYER_SHOT_TAG, MPI_COMM_WORLD, &flag, &status);
            if (flag) {
                int shot;
                MPI_Recv(&shot, 1, MPI_INT, MASTER, PLAYER_SHOT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                int random = rand() % 100;
                if (random < 20){
                    //shoot left
                    int kill = 1;
                    MPI_Send(&kill, 1, MPI_INT, left, KILL_TAG, MPI_COMM_WORLD);
                } else if (random < 35){
                    //shoot right
                    int kill = 1;
                    MPI_Send(&kill, 1, MPI_INT, right, KILL_TAG, MPI_COMM_WORLD);
                }else if (random < 55){
                    //block
                    alive = 1;
                }else{
                    //get hit
                    alive = 0;
                }
            }

            //Invader has been killed
            if (left != MPI_PROC_NULL){
                MPI_Iprobe(left, KILL_TAG, MPI_COMM_WORLD, &flag, &status);
                if (flag) {
                    int kill;
                    MPI_Recv(&kill, 1, MPI_INT, left, KILL_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    alive = 0;
                }
            }

            if (right != MPI_PROC_NULL){
                MPI_Iprobe(right, KILL_TAG, MPI_COMM_WORLD, &flag, &status);
                if (flag) {
                    int kill;
                    MPI_Recv(&kill, 1, MPI_INT, right, KILL_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    alive = 0;
                }
            }
            // 4. Check for game over signal
            MPI_Iprobe(MASTER, 99, MPI_COMM_WORLD, &flag, &status);
            if (flag) {
                int stop;
                MPI_Recv(&stop, 1, MPI_INT, MASTER, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                break;
            }
            
            // 5. Synchronize with master
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    
    MPI_Finalize();
    return 0;
}

void print_game(int n, int m, Invader* invaders, int player_col, 
                Cannonball* player_balls, int num_pballs,
                Cannonball* enemy_balls, int num_eballs) {
    printf("\n--- Game Board ---\n");
    

    for (int row = 0; row < n; row++) {
        for (int col = 0; col < m; col++) {
            int index = row * m + col;
            if (invaders[index].alive) {
                printf("[👾 ]");
            } else {
                printf("[ ]");
            }
            printf(" ");
        }
        printf("\n");
    }
    

    printf("\n");
    
    
    for (int col = 0; col < m; col++) {
        if (col == player_col) {
            printf(" ^  ");
        } else {
            printf("    ");
        }
    }
    printf(" <- Player\n");
    printf("------------------\n");
}
