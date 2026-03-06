from collections import defaultdict
import time
import argparse

def getBCumiDict(rawEvent):
    BCumi_dict = defaultdict(lambda: defaultdict(list))
    for each in rawEvent.split('>')[1:]:
        name = each.split('\n')[0]
        spl = name.strip().split(" ")
        BC = spl[1]
        UMI = spl[2]
        length = spl[4].split("=")[1]
        numIndel = spl[11].split("=")[1]
        edit_mutation = each.split('\n')[1]
        BCumi_dict[BC][UMI] = [numIndel, length, name, edit_mutation]
    return BCumi_dict
    
def getMinIndel(BCumiInfo_dict):
    result_dict = {}
    for B_, u_info in BCumiInfo_dict.items():
        result = []
        for i in u_info.keys():
            #if int(u_info[i][0]) <= 5:
            if int(u_info[i][0]) <= 20:
            #if int(u_info[i][0]) <= 100:
                result.append(f'>{u_info[i][2]}\n{u_info[i][3]}')

        result_dict[B_] = result
    return result_dict
    
    
def Parsers():
    parser = argparse.ArgumentParser(description="Get Events base on minIndels")
    parser.add_argument('-i', '--input', type=str, nargs='?' ,required=True, help="raw Event file")
    parser.add_argument('-o', '--outfile', type=str, nargs='?', required=True, help="Output editevents file")
    args = parser.parse_args()
    return args    
       
if __name__ == "__main__":
    start = time.time()
    options = Parsers()
    
    inRawEvent = open(options.input,'r').read()
    BC_umi_info_dict = getBCumiDict(inRawEvent)
    results_dict = getMinIndel(BC_umi_info_dict)
    outEvent = open(options.outfile, 'w')
    for key, value in results_dict.items():
        for item in value:
            outEvent.write(item + '\n')
    outEvent.close()
     
    end = time.time()
    print(end - start)
    
    