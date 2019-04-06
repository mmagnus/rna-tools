def seq_numbers():
    """
    output:
    - |0.......|10......|20......|30......|40......|50......|60......|70......|80......|90......|0.......|10......|20......|30......|40......|50......|60......|70
    """
    c=0
    f=c*100
    s=c*100+100
    lineNr='|1.......'#beginning of the line
    for i in [10,20,30,40,50,60,70,80,90]:
            lineNr+= '|' + str(i)+('.'*(9-len(str(i))))
    lineNr = lineNr + lineNr + lineNr
    return lineNr

if __name__ == '__main__':
    print((seq_numbers()))
    
