from itertools import combinations_with_replacement
sequences = ['ATTCG', 'ATTAACG', 'ATGCG']


dico_score = {
    'AA': 3,
    'CC': 3,
    'GG': 3,
    'TT': 3,
    'AG': 1,
    'CT': 1,
    'AT': 0,
    'AC': 0,
    'CG': 0,
    'GT': 0,
    '-A': 0,
    '-C': 0,
    '-G': 0,
    '-T': 0,
    '--': 0,
}

#Fonction permettant d'insérer les gaps pour donner la séquence correspondant à l'alignement devant être testé
def insert_dash(chaine, pos):
  return chaine[:pos]+'-'+chaine[pos:]

#Fonction de création d'un alignement à partir d'une séquence et de positions de gap
def create_alignment(seq, gap):
  new_seq = seq
  for i in range(len(gap)):
    new_seq = insert_dash(new_seq, gap[i]+i)
  return new_seq

#Fonction calculant le score d'alignement d'un alignement multiple de 3 séquences
def score_calcul(alignment):
  gap_num = 0
  gap_num += alignment[0].count('-')
  gap_num += alignment[1].count('-')
  gap_num += alignment[2].count('-')
  score = 0
  for i in range(len(alignment[0])):
    if alignment[0][i] > alignment[1][i]:
      score += dico_score[alignment[1][i] + alignment[0][i]]
    else:
      score += dico_score[alignment[0][i] + alignment[1][i]]
    if alignment[0][i] > alignment[2][i]:
      score += dico_score[alignment[2][i] + alignment[0][i]]
    else:
      score += dico_score[alignment[0][i] + alignment[2][i]]
    if alignment[1][i] > alignment[2][i]:
      score += dico_score[alignment[2][i] + alignment[1][i]]
    else:
      score += dico_score[alignment[1][i] + alignment[2][i]]
  score -= gap_num*3
  return score
score_max = 0
align_max = ''
min_len = len(sequences[0])

#Estimation du nombre maximale de gap insérable dans l'alignement
for seq in sequences:
  if len(seq) > min_len:
    min_len = len(seq)
max_len = int(min_len*1.5) + (min_len % 2) - 1


counter = 0
#Calcul de tous les arrangements de gap possibles dans les alignements
for ngaps in range(min_len,max_len+1):
  print(f"nagps : {ngaps}")
  for i in combinations_with_replacement(range(len(sequences[0])+1),ngaps - len(sequences[0])):
    for j in combinations_with_replacement(range(len(sequences[2])+1),ngaps - len(sequences[2])):
      for k in combinations_with_replacement(range(len(sequences[1])+1),ngaps - len(sequences[1])):
        counter += 1
        sequence_gen = [create_alignment(sequences[0], i), create_alignment(sequences[1], k), create_alignment(sequences[2], j)]
        score = score_calcul(sequence_gen)
        if counter % 100000 == 0:
          print(f'compte: {counter}, score: {score}')
        if score > score_max:
          score_max = score
          align_max = sequence_gen
          print(f'score : {score_max}')
          print(f'{align_max[0]}\n{align_max[1]}\n{align_max[2]}')
print(f'compte :{counter}')
