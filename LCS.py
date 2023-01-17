import model
import pickle

matrix = {
    'indel': 0.01,
    'mut': 0.1
}

list_results = []

#Double boucle sur la taille et le nombre des séquences

for i in range(2,21):
  #Boucle sur la taille des séquences
  for j in range (2,21):
    #Boucle sur le nombre de séquence
    #Premier duplicat
    print(i,j,1)
    #Génération des séquences étudiées
    seqs = model.seq_generator_mutate(alphabet = 'CSTAGPDEQNHRKMILVWYF', seq_size = i, seq_number=j, muta_matrix=matrix)
    LCS = model.longest_common_subsequence(seqs)
    LCS['seqs']=seqs
    LCS['nb_seq'] = j
    LCS['size_seq'] = i
    list_results.append(LCS)
    #Second duplicat
    print(i,j,2)
    seqs = model.seq_generator_mutate(alphabet = 'CSTAGPDEQNHRKMILVWYF', seq_size = i, seq_number=j, muta_matrix=matrix)
    LCS = model.longest_common_subsequence(seqs)
    LCS['seqs']=seqs
    LCS['nb_seq'] = j
    LCS['size_seq'] = i
    list_results.append(LCS)
  with open("result_20_20", 'wb') as f1:
    pickle.dump(list_results, f1)