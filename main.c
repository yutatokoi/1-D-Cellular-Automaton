/* Program to perform 1D cellular automaton (CA) computations and to use 1D CA
   to solve the density classification problem.

  Skeleton program written by Artem Polyvyanyy, http://polyvyanyy.com/,
  September 2024, with the intention that it be modified by students
  to add functionality, as required by the assignment specification.
  All included code is (c) Copyright University of Melbourne, 2024.
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

/* #DEFINE'S -----------------------------------------------------------------*/
#define SDELIM "==STAGE %d============================\n"   // stage delimiter
#define MDELIM "-------------------------------------\n"    // delimiter of -'s
#define THEEND "==THE END============================\n"    // end message
#define POSSIBLE_NEIGHBOURHOODS " 000 001 010 011 100 101 110 111\n"

#define CRTRNC '\r'      // carriage return character
#define NBRHDS 8         // number of possible neighborhoods
#define CELLS_IN_NBRHD 3 // number of cells in a neighborhood
#define ON_STATE '*'     // char representing the on state of a cell
#define OFF_STATE '.'    // char representing the off state of a cell
#define OFF_INT '0'      // char representing the off state of a cell as a num

#define RULE_0 "00000000" // rule_t representing update rule 0

#define RULE184 184
#define RULE232 232

#define LESS_THAN '<'
#define MORE_THAN '>'
#define EQUAL     '='

/* TYPE DEFINITIONS ----------------------------------------------------------*/
typedef char cells_t;                   // base type to store states of cells
typedef struct state state_t;           // a cellular automaton state
typedef unsigned char rule_t[NBRHDS+1]; // a CA update rule

struct state {                   // a state in a CA is defined by
    cells_t*        clls;        // ... an array of cells and
    state_t*        next;        // ... a link to the next state
};

typedef struct {                 // a run of a CA consists of
    state_t*        init;        // ... the initial state and
    state_t*        curr;        // ... the current state,
} run_t;                         // implemented as a linked list of states

typedef struct {                 // an elementary CA is defined by
    unsigned int    size;        // ... a number of cells,
    unsigned int    time;        // ... the current time step,
    rule_t          rule;        // ... an update rule function, and
    run_t*          run;         // ... a run of state steps
    unsigned int    rule_int;    // ... int representation of the update rule
} CA_t;

/* FUNCTION PROTOTYPES -------------------------------------------------------*/

/* USEFUL FUNCTIONS ----------------------------------------------------------*/
int mygetchar(void);             // getchar() that skips carriage returns
void int_to_rule(int input, rule_t *p);

/* Cellular automaton functions ----------------------------------------------*/
void print_cells(cells_t *cells, int time);
void print_generations(CA_t ca, int start_time, int steps);
void execute_automaton(CA_t ca, int rule, int start_time, int steps);
int count_on_states(CA_t ca, int cell, int start_time);
void free_ca(CA_t* ca);

/* Stage 0 Functions ---------------------------------------------------------*/
CA_t* read_ca();
void print_stage0(CA_t ca);

/* Stage 1 Functions ---------------------------------------------------------*/
int nbrhd_to_int(cells_t *nbrhd);
cells_t int_to_cell(int num);
cells_t* next_generation(cells_t *state, int size, rule_t rule);
void print_stage1(CA_t ca, int stage1_time, int stage1_cell, int stage1_start);

/* Stage 2 Functions ---------------------------------------------------------*/
char final_state(CA_t ca);
void print_stage2(CA_t ca, int stage1_time, int stage2_cell, int stage2_start,
                  int n_steps, int m_steps);

/* WHERE IT ALL HAPPENS ------------------------------------------------------*/
int main(int argc, char *argv[]) {
    int stage = 0;
    printf(SDELIM, stage++);
    CA_t* ca = read_ca(); 
    
    // Read from stdin lines 4~7 of the prescribed input
    int stage1_time;
    scanf("%d\n", &stage1_time);
    assert(0 <= stage1_time);
    int stage1_cell, stage1_start;
    scanf("%d,%d\n", &stage1_cell, &stage1_start);
    assert((0 <= stage1_cell) && (stage1_cell <= ca->size));
    assert((0 <= stage1_start) && (stage1_start <= stage1_time));
    int stage2_cell, stage2_start;
    scanf("%d,%d\n", &stage2_cell, &stage2_start);
    assert((0 <= stage2_cell) && (stage2_cell <= ca->size));
    int n_steps = (int)floor((ca->size - 2) / 2),
        m_steps = (int)floor((ca->size - 1) / 2);
    assert((0 <= stage2_start) && (stage2_start <= (n_steps + m_steps + 
                                   stage1_time)));

    print_stage0(*ca);

    printf(SDELIM, stage++);
    if (stage1_time > 0) {
        execute_automaton(*ca, ca->rule_int, 0, stage1_time);
    }
    print_stage1(*ca, stage1_time, stage1_cell, stage1_start);

    printf(SDELIM, stage);
    execute_automaton(*ca, RULE184, stage1_time, n_steps);
    execute_automaton(*ca, RULE232, (stage1_time + n_steps), m_steps);
    print_stage2(*ca, stage1_time, stage2_cell, stage2_start, n_steps, m_steps);

    free_ca(ca);
    printf(THEEND);
    return EXIT_SUCCESS;
}

/* USEFUL FUNCTIONS ----------------------------------------------------------*/

// An improved version of getchar(); skips carriage return characters.
// NB: Adapted version of the mygetchar() function by Alistair Moffat
int mygetchar() {
    int c;
    while ((c=getchar())==CRTRNC);          // skip carriage return characters
    return c;
}

/* Convert an integer representation of a rule code (0~255) to its binary 
   representation as a string. Ex: rule code 30 converts to:
   
   111 110 101 100 011 010 001 000 
    0   0   0   1   1   1   1   0
   -> "00011110" 
*/
void int_to_rule(int input, rule_t *p) {
    rule_t rule = RULE_0;
    for (int i = NBRHDS - 1; i >= 0; i--) {
        if (input - (int)pow(2.0, (double)i) >= 0) {
            input -= (int)pow(2.0, (double)i);
            rule[NBRHDS-i-1] = '1';
        }
    }
    strcpy((char *)*p, (char *)rule);
}

/* Cellular automaton functions ----------------------------------------------*/

void print_cells(cells_t *cells, int time) {
    printf("%4d: %s\n", time, cells);
}

void print_generations(CA_t ca, int start_time, int steps) {
    state_t *curr_state = ca.run->init;
    int t = 0;
    // Traverse linked list to get to the starting time
    while (t < start_time) {
        curr_state = curr_state->next;
        t++;
    }

    // Print out the clls starting from time t, until `steps` number of steps
    // have been traversed
    while (t <= start_time + steps) {
        printf("%4d: %s\n", t, curr_state->clls);
        curr_state = curr_state->next;
        t++;
    }
}

// Execute the evolution algorithm as specified by the generation rule
void execute_automaton(CA_t ca, int rule, int start_time, int steps) {
    state_t *curr_state = ca.run->init;
    // Traverse linked list to get to the starting time
    while (start_time > 0) {
        curr_state = curr_state->next;
        start_time--;
    }

    rule_t rule_code;
    int_to_rule(rule, &rule_code);

    // Generate and store the next generation for `steps` number of times, using
    // the rule number `rule`
    while (steps > 0) {
        state_t *next_state = malloc(sizeof(state_t));
        assert(next_state != NULL);
        next_state->clls = next_generation(curr_state->clls, ca.size, 
                                           rule_code);
        next_state->next = NULL;
        curr_state->next = next_state;
        curr_state = next_state;
        steps--;
    }
    ca.run->curr = curr_state;
}

// Return the number of times the specified cell is in the on state, starting
// from t = `start_time`, until the end of the linked list of states.
int count_on_states(CA_t ca, int cell, int start_time) {
    int on = 0;
    state_t *curr_state = ca.run->init;
    
    // Ignore the lines until the clls at t = start_time
    while (start_time > 0) {
        curr_state = curr_state->next;    
        start_time--;
    }

    // Count number of times the cell is in the on state until the end of the 
    // run is reached.
    while (curr_state != NULL) {
        if (curr_state->clls[cell] == ON_STATE) on++;
        curr_state = curr_state->next;
    }

    return on;
}

// Function to free all malloc'd memory at the end of execution
void free_ca(CA_t* ca) {
    state_t *current_state = ca->run->init;
    while (current_state != NULL) {
        state_t *next_state = current_state->next;
        free(current_state->clls);
        free(current_state);
        current_state = next_state;
    }
    free(ca->run);
    free(ca);
}

/* Stage 0 Functions ---------------------------------------------------------*/

// Read the first 3 lines from stdin that describes the initial state and rule
// of the cellular automaton
CA_t* read_ca() {
    CA_t* ca = malloc(sizeof(CA_t));
    assert(ca != NULL);
    ca->time = 0;
    scanf("%d\n", &ca->size);
    scanf("%d\n", &ca->rule_int);
    int_to_rule(ca->rule_int, &ca->rule);

    state_t *state = malloc(sizeof(state_t));
    assert(state != NULL);

    // Read the string representing the initial state of the cellular automaton
    size_t size = ca->size;
    cells_t *cells = malloc((size + 1) * sizeof(cells_t));
    assert(cells != NULL);

    int ch;
    size_t length = 0;
    while ((ch = mygetchar()) != '\n') {
        assert(length < size);
        cells[length++] = (cells_t)ch;
    }
    cells[length] = '\0';

    state->clls = cells;
    state->next = NULL;

    run_t *run = malloc(sizeof(run_t));
    assert(run != NULL);
    run->init = state;
    run->curr = state;

    ca->run = run;

    return ca;
}

void print_stage0(CA_t ca) {
    printf("SIZE: %d\n", ca.size);
    printf("RULE: %d\n", ca.rule_int);
    printf(MDELIM);
    printf(POSSIBLE_NEIGHBOURHOODS);
    for (int i = NBRHDS - 1; i >= 0; i--) {
        if (i == NBRHDS - 1) {
            printf("%3c", ca.rule[i]);
        } else if (i == 0) {
            printf("%4c \n", ca.rule[i]);
        } else {
            printf("%4c", ca.rule[i]);
        }
    }
    printf(MDELIM);
    print_cells(ca.run->init->clls, 0);
}

/* Stage 1 Functions ---------------------------------------------------------*/

// Convert a cells_t* representing a neighborhood, into the corresponding
// decimal value. Ex: "*.*" is 0b101, so will return 5 in decimal
int nbrhd_to_int(cells_t *nbrhd) {
    int rule_int = 0;
    int exponent = CELLS_IN_NBRHD - 1;
    for (int i = 0; i < CELLS_IN_NBRHD; i++) {
        int cell_state = (nbrhd[i] == ON_STATE);
        rule_int += cell_state * (int)pow(2.0, exponent--);
    }
    return rule_int;
}

// Convert 0 (off) to '.', 1 (on) to '*'
cells_t int_to_cell(int num) {
    if (num == 0) {
        return OFF_STATE;
    } else if (num == 1) {
        return ON_STATE;
    } else {
        exit(EXIT_FAILURE);
    }
}

// Given a state of the cellular automaton, return the next generation
cells_t* next_generation(cells_t *state, int size, rule_t rule) {
    cells_t *new_state = malloc((size + 1) * sizeof(cells_t));
    assert(new_state != NULL);
    for (int i = 0; i < size; i++) {
        int j = i - 1;
        int nbrhd_i = 0;
        cells_t nbrhd[CELLS_IN_NBRHD+1];
        cells_t new_cell;
        if (i == 0) {
            // Wrap to the end
            nbrhd[nbrhd_i++] = state[size-1];
            nbrhd[nbrhd_i++] = state[++j];
            nbrhd[nbrhd_i] = state[++j];
        } else if (i == size - 1) {
            // Wrap to the beginning 
            nbrhd[nbrhd_i++] = state[j++];
            nbrhd[nbrhd_i++] = state[j];
            nbrhd[nbrhd_i] = state[0];
        } else {
            nbrhd[nbrhd_i++] = state[j++];
            nbrhd[nbrhd_i++] = state[j++];
            nbrhd[nbrhd_i] = state[j];
        }
        new_cell = int_to_cell(rule[NBRHDS - 1 - nbrhd_to_int(nbrhd)]
                               - OFF_INT);
        new_state[i] = new_cell;
    }
    new_state[size] = '\0';
    return new_state;
}

void print_stage1(CA_t ca, int stage1_time, int stage1_cell, int stage1_start) {
    state_t *curr = ca.run->init;
    // Print all the states that were generated in Stage 1
    for (int t = 0; t <= stage1_time; t++) {
        print_cells(curr->clls, t);
        curr = curr->next;
    }
    printf(MDELIM);

    int on = count_on_states(ca, stage1_cell, stage1_start);
    int off = stage1_time - stage1_start - on + 1;
    printf("#ON=%d #OFF=%d CELL#%d START@%d\n", on, off, 
                                                stage1_cell, stage1_start);
}

/* Stage 2 Functions ---------------------------------------------------------*/

// If all the cells in the final state are ON, return '>'
// If all the cells in the final state are OFF, return '<'
// If there are equal number of ON and OFF cells, return '='
char final_state(CA_t ca) {
    state_t *state = ca.run->init;

    // Traverse the linked list until hitting the final state
    while (state->next != NULL) {
        state = state->next;
    }

    cells_t *clls = state->clls;
    if ((clls[0] == ON_STATE) && (clls[ca.size - 1] == ON_STATE)) {
        return MORE_THAN;
    } else if ((clls[0] == OFF_STATE) && (clls[ca.size - 1] == OFF_STATE)) {
        return LESS_THAN;
    } else {
        return EQUAL;
    }
}

void print_stage2(CA_t ca, int stage1_time, int stage2_cell, int stage2_start,
                  int n_steps, int m_steps) {
    printf("RULE: %d; STEPS: %d.\n", RULE184, n_steps);
    printf(MDELIM);
    print_generations(ca, stage1_time, n_steps);
    printf(MDELIM);
    printf("RULE: %d; STEPS: %d.\n", RULE232, m_steps);
    printf(MDELIM);
    print_generations(ca, (stage1_time + n_steps), m_steps);
    printf(MDELIM);
    int on = count_on_states(ca, stage2_cell, stage2_start);
    int off = stage1_time + n_steps + m_steps - stage2_start - on + 1;
    printf("#ON=%d #OFF=%d CELL#%d START@%d\n", on, off, 
                                                stage2_cell, stage2_start);
    printf(MDELIM);
    
    state_t *curr_state = ca.run->init;
    int t = 0;
    while (t < stage1_time) {
        curr_state = curr_state->next;
        t++;
    }
    print_cells(curr_state->clls, t);

    printf("AT T=%d: #ON/#CELLS %c 1/2\n", stage1_time, final_state(ca));
}
