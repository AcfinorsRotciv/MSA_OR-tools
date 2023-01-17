from ortools.linear_solver import pywraplp
import numpy as np
import copy
import pickle

def seq_generator(alphabet, seq_size, seq_number):
  gen_list = []
  for i in range(seq_number):
    seq = ''.join(list(np.random.choice(list(alphabet), seq_size)))
    gen_list.append(seq)
  return gen_list

def seq_generator_mutate(alphabet, seq_size, seq_number, muta_matrix):
  gen_list = []
  taille = seq_size
  seq = ''.join(list(np.random.choice(list(alphabet), taille)))
  for i in range(seq_number):
    seq_new = []
    for pos in range(len(seq)):
      rando = np.random.rand()
      if rando < muta_matrix['indel']:
        if np.random.rand() < 0.5:
          seq_new.append(np.random.choice(list(alphabet)))
          pos -= 1
      else:
        if rando < muta_matrix['mut']:
          seq_new.append(np.random.choice(list(alphabet)))
        else:
          seq_new.append(seq[pos])
    gen_list.append(''.join(seq_new))
  return gen_list

def longest_common_subsequence(sequences):
    """
    Finds the longest common subsequence of a list of sequences.

    Args:
        sequences: A list of sequences.

    Returns:
        A list representing the longest common subsequence of the sequences.
    """
    # Solver
    # Create the mip solver with the SCIP backend.
    solver = pywraplp.Solver.CreateSolver('SCIP')

    if not solver:
        return

    # Create a solver
    # solver = pywraplp.Solver('SolveIntegerProblem',
    #                          pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)

    # Define the decision variables and reading the data

    # declaration of match variable m and varaible encoding the fact that the match belongs to the LCS
    # mi,j,p,q is set to 1 for each pair of consecutive sequences si and sj if si,p = sj,p 0 else.
    m = []
    x = {}
    for i in range(len(sequences)):
      if i < (len(sequences)-1):
        m_ij = []
        for p in range(len(sequences[i])):
          m_ijp = []
          for q in range(len(sequences[i+1])):
            m_ijp.append(1) if sequences[i][p] == sequences[i+1][q] else m_ijp.append(0)
            x[i, p, q] = solver.IntVar(0, 1, f'x_i{i}_j{i+1}_p{p}_q{q}')
          m_ij.append(m_ijp)
        m.append(m_ij)

    # Define the constraints

    # Constraint #1

    # Symbols si[p] and sj[q] can be part of a longest common 
    # subsequence (xi,j,p,q =1) if and only if they are identical
    for i in range(len(m)):
      for p in range(len(m[i])):
        for q in range(len(m[i][p])):
          solver.Add(x[i, p, q] <= m[i][p][q])

    # Constraint #2

    # A symbol si[p] can be aligned with at most one symbol sj[q]
    for i in range(len(m)):
      for p in range(len(m[i][0])):
        solver.Add(solver.Sum([(x[i, q, p]) for q in range(len(m[i]))]) <= 1)
    
    # Constraint #3

    # A symbol sj[q] can be aligned with at most one symbol si[p]
    for i in range(len(m)):
      for p in range(len(m[i])):
        solver.Add(solver.Sum([(x[i, p, q]) for q in range(len(m[i][p]))]) <= 1)

    # Constraint #4

    # A symbol si[p] that is aligned
    # with a symbol sj [q] must also
    # be aligned with some symbol sk[r].
    # For each three consecutive sequences
    # si, sj and sk.
    # |sk|
    #  ∑ xj,k,q,r >= xi,j,p,q
    # r=1
    for i in range(len(m)):
      if i < (len(m)-1):
        j = i+1
        for p in range(len(m[i])):
          for q in range(len(m[i][p])):
            solver.Add(solver.Sum([x[j, q, r] for r in range(len(m[j][q]))]) >= x[i, p, q])

    # A symbol si[p] that is aligned
    # with a symbol sj [q] must also
    # be aligned with some symbol sk[r].
    # For each three consecutive sequences
    # si, sj and sk.
    # |si|
    #  ∑ xi,j,p,q >= xj,k,q,r
    #  p=1
    for i in range(len(m)):
      if i < (len(m)-1):
        for q in range(len(m[i+1])):
          for r in range(len(m[i+1][q])):
            solver.Add(solver.Sum([x[i, p, q] for p in range(len(m[i]))]) >= x[i+1, q, r])

    # Constraint #5

    # The symbols that comprise a solution
    # must appear in the same order in all the
    # sequences. This means, there cannot be
    # “crosses” between these symbols in any
    # pair of consecutive sequences.
    for i in range(len(m)):
      for p1 in range(len(m[i])):
        for p2 in range(len(m[i])):
          for q1 in range(len(m[i][p1])):
            for q2 in range(len(m[i][p1])):
              if (((p1 > p2) and (q1 < q2)) or ((p1 < p2) and (q1 > q2))):
                solver.Add(solver.Sum([x[i, p1, q1], x[i, p2, q2]]) <= 1)


    # Define the objective function
    # The objective is to maximize the length of the longest common subsequence
    solver.Maximize(solver.Sum([x[i, p, q] for i in range(len(m)) for p in range(len(m[i])) for q in range(len(m[i][p]))]))
    
    # Solve the mixed integer program
    status = solver.Solve()

    # Check the status of the solution
    if status == pywraplp.Solver.OPTIMAL:

      print('\nAdvanced usage:')
      print('Problem solved in %f seconds' % (solver.wall_time()/1000) )
      print('Problem solved in %d iterations' % solver.iterations())
      print('Problem solved in %d branch-and-bound nodes' % solver.nodes())

        
      LCS = ""
      for i in range(len(m)):
        for p in range(len(m[i])):
          for q in range(len(m[i][p])):
            if x[i, p, q].solution_value() == 1:
              # print(f'{x[i, p, q].name()} = {x[i, p, q].solution_value()}')
              if i == 0: LCS += sequences[i][p]
            # test si les solutions respectent les contraintes
            # if m[i][p][q] == 1:
            #   print(f'm_i{i}_j{i+1}_p{p}_q{q} = {m[i][p][q]}')
      retour = {}
      retour['time']=solver.wall_time()/1000
      retour['iter']=solver.iterations()
      retour['nodes']=solver.nodes()
      retour['seq']=LCS
    else:
      print('The problem does not have an optimal solution.')

    return retour

with open("BLOSUM62", "r") as f:
  index = f.readline().strip("\n").replace("  ", " ").strip().split(" ")
  BLOSUM = {}
  for line in f:
    line_list = (line.strip("\n").replace("  ", " ").strip().split(" "))
    if len(line_list) > 1:
      for i, AA in enumerate(index):
        BLOSUM[AA+line_list[0]] = int(line_list[i+1])

def multiple_sequence_alignment_ILP(sequences, gap_penalty):
    """
    Finds the longest common subsequence of a list of sequences.

    Args:
        sequences: A list of sequences.

    Returns:
        A list representing the longest common subsequence of the sequences.
    """
    # Solver
    # Create the mip solver with the SCIP backend.
    solver = pywraplp.Solver.CreateSolver('SCIP')
    if not solver:
        return

    gap_penalty = -4
    # Solver
    solver = pywraplp.Solver('MultipleSequenceAlignment',
                              pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)
    
    solver.SetTimeLimit(600*1000)
    # Define the decision variables
    x = {}
    # binary variables, x_i,j,p,q = 1 if symbols si[p] and sj[q] are part of a multiple alignment
    for i in range(len(sequences) - 1):
      j = i+1
      for p in range(len(sequences[i])):
        for q in range(len(sequences[j])):
          x[i, j, p, q] = solver.IntVar(0, 1, f'x_{i}_{j}_{p}_{q}')


    # Define the constraints
    # Constraint 1: Symbols si[p] and sj[q] can be part of a longest common subsequence (xi,j,p,q = 1) if and only if they are identical (mi,j,p,q = 1). 
    # for i in range(len(sequences)):
    #     for j in range(i+1,len(sequences)):
    #         for p in range(len(sequences[i])):
    #             for q in range(len(sequences[j])):
    #               solver.Add(x[i, j, p, q] <= blosum_matrix[sequences[i][p]][sequences[j][q]])

    # Constraint 2: A symbol si[p] can be aligned with at most one symbol sj[q].
    for i in range(len(sequences) - 1):
      j = i+1
      for q in range(len(sequences[j])):
        solver.Add(solver.Sum([x[i, j, p, q] for p in range(len(sequences[i]))]) <= 1)

    # Constraint 3: A symbol sj[q] can be aligned with at most one symbol si[p].
    for i in range(len(sequences) - 1):
      j = i+1
      for p in range(len(sequences[i])):
        solver.Add(solver.Sum([x[i, j, p, q] for q in range(len(sequences[j]))]) <= 1)
                
    # Constraint 4: A symbol si[p] that is aligned with a symbol sj[q] must also be aligned with some symbol sk[r].
    for i in range(len(sequences) - 2):
      j = i+1
      k = i+2
      for p in range(len(sequences[i])):
        for q in range(len(sequences[j])):
          solver.Add(x[i, j, p, q] <= solver.Sum([x[j, k, q, r] for r in range(len(sequences[k]))]))

    for i in range(len(sequences) - 2):
      j = i+1
      k = i+2
      for q in range(len(sequences[j])):
        for r in range(len(sequences[k])):
          solver.Add(x[j, k, q, r] <= solver.Sum([x[i, j, p, q] for p in range(len(sequences[i]))]))

    # Constraint 5: The symbols that comprise a solution must appear in the same order in all the sequences.
    for i in range(len(sequences) - 1):
      j = i+1
      for p1 in range(len(sequences[i])):
        for q1 in range(len(sequences[j])):
          for p2 in range(len(sequences[i])):
            for q2 in range(len(sequences[j])):
              if (((p1 > p2) and (q1 < q2)) or ((p1 < p2) and (q1 > q2))):
                solver.Add(solver.Sum([x[i,j, p1, q1], x[i,j, p2, q2]]) <= 1)


    # Define the objective function
    # Maximize the alignment score

    objective = [solver.Sum([x[i, i+1, p, q] * BLOSUM[sequences[i][p] + sequences[i+1][q]] for i in range(len(sequences) - 1) for p in range(len(sequences[i])) for q in range(len(sequences[i+1]))])]
    objective.append(gap_penalty * solver.Sum([1 - solver.Sum([x[i, i+1, p, q] for q in range(len(sequences[i+1]))]) for i in range(len(sequences) - 1) for p in range(len(sequences[i]))]))
    objective.append(gap_penalty * solver.Sum([1 - solver.Sum([x[i, i+1, p, q] for p in range(len(sequences[i]))]) for i in range(len(sequences) - 1) for q in range(len(sequences[i+1]))]))
    solver.Maximize(solver.Sum(objective))

    # Solve the problem
    status = solver.Solve()

    if status == pywraplp.Solver.OPTIMAL:
      # Extract the final alignment from the solution
      print('Problem solved in %f seconds' % (solver.wall_time()/1000) )
      print('Problem solved in %d iterations' % solver.iterations())
      print('Problem solved in %d branch-and-bound nodes' % solver.nodes())
      print(solver.Objective().Value())
      
      max_seq = 0
      for i in sequences:
        if (len(i) > max_seq):
          max_seq = len(i) 

      MSA = []
      for p in range(max_seq):
        for i in range(len(sequences) - 1):
          j = i+1
          no_match = True
          for q in range(max_seq):
            if (p < len(sequences[i])) and (q < len(sequences[j])):
              if x[i, j, p, q].solution_value() == 1:
                MSA.append((f"{sequences[i][p]+sequences[j][q]}", (p, q), (i, j)))
                no_match = not no_match
              elif ((q == (len(sequences[j]) - 2)) and (x[i, j, p, q+1].solution_value() == 0)) and no_match:
                MSA.append((f"{sequences[i][p]}-", (p, None), (i, j)))

      
      retour = {}
      retour['time']=solver.wall_time()/1000
      retour['iter']=solver.iterations()
      retour['nodes']=solver.nodes()
      retour['align']=MSA
      return retour
      
    else:
        print("No solution found")
        retour = {}
        retour['time']=solver.wall_time()/1000
        retour['iter']=solver.iterations()
        retour['nodes']=solver.nodes()
        retour['align']='NA'
        return retour
