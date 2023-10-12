import json
import pandas as pd
import os
from collections import defaultdict
from tqdm import tqdm


def process_gene():
    # gene_symbol_t = defaultdict(list)
    # gene_symbol = {'P16112': ['ACAN'], 'O75564': ['JRK'], 'Q9UFH2': ['DNAH17'], 'P03999': ['OPN1SW'], 'P16455': ['MGMT'],
    #                'Q96MI9': ['AGBL1'], 'Q68DE3': ['USF3'], 'P19013': ['KRT4'], 'Q8IWN7': ['RP1L1'], 'Q70YC4': ['ZNF365'],
    #                'P0DPH7': ['TUBA3C'], 'P0DPH8': ['TUBA3D']}
    symbol = []
    with open('../data/disgenet/gene_symbol_t.json', 'r') as f:
        data = json.load(f)
    # with open('gene_uniprot.json', 'r') as fu:
    #     gene_un = json.load(fu)
    # print(len(gene_un))
    # print(len(list(data.keys())))
    # for j in list(data.keys()):
    #     for i, k in gene_un.items():
    #         if j in k:
    #             gene_symbol_t[j].append(i)
    # for z, v in gene_symbol.items():
    #     gene_symbol_t[z] = gene_symbol_t[z] + v
    # print(gene_symbol_t)
    # print(len(gene_symbol_t))
    # num = 0
    for z, v in data.items():
        symbol = symbol + v
    print(len(symbol))
    with open('symbol.json', 'w') as fg:
        json.dump(list(set(symbol)), fg, indent=2)


def process_disease_hpo():
    data = pd.read_csv('../data/disgenet/disease_mappings.tsv', sep='\t')
    print(data[['diseaseId']].values.tolist())
    DisGeNet_disease_hpo = defaultdict(list)
    number = 0
    for d, t, h in data[['diseaseId', 'vocabulary', 'code']].values.tolist():
        if t == 'HPO':
            number = number + 1
            DisGeNet_disease_hpo[d].append(h)
    print(number)
    num = 0
    for j, k in DisGeNet_disease_hpo.items():
        DisGeNet_disease_hpo[j] = list(set(k))
        num = num + len(DisGeNet_disease_hpo[j])
    print(len(DisGeNet_disease_hpo))
    print(num)
    # with open('DisGeNet_disease_hpo.json', 'w') as fg:
    #     json.dump(DisGeNet_disease_hpo, fg, indent=2)


def process_gene_hpo():
    symbol_hpo = defaultdict(list)
    no_list = []
    with open('../data/disgenet/symbol.json', 'r') as f:
        data = json.load(f)
    print(len(data))
    with open('../data/disgenet/DisGeNet_disease_hpo.json', 'r') as fg:
        DisGeNet_disease_hpo = json.load(fg)
    for i in data:
        path = 'gene_disease/{}.csv'.format(i)
        if os.path.exists(path):
            data1 = pd.read_csv(path, sep='\t')

            for j in data1.values:
                # print(j)
                key = j[0].replace('"', '').split(',')[1]
                # print(key)
                if key in DisGeNet_disease_hpo.keys():
                    symbol_hpo[i] = list(set(symbol_hpo[i] + DisGeNet_disease_hpo[key]))
        else:
            no_list.append(i)
    with open('../data/disgenet/no_list.json', 'w') as fgg:
        json.dump(no_list, fgg, indent=2)
    print(len(symbol_hpo))
    num = 0
    zero_num = []
    for j, i in symbol_hpo.items():
        if len(i) > 0:
            num = num + len(i)
        else:
            zero_num.append(j)
    print(num)
    print(len(zero_num))
    with open('../data/disgenet/DisGeNet_gene_hpo.json', 'w') as fggg:
        json.dump(symbol_hpo, fggg, indent=2)


def un_hpo():
    un_hpo_t = defaultdict(list)
    with open('../data/disgenet/DisGeNet_gene_hpo.json', 'r') as f:
        data = json.load(f)
    with open('../data/disgenet/gene_symbol_t.json', 'r') as fg:
        gene_symbol_t = json.load(fg)
    print(len(gene_symbol_t))
    for i, j in gene_symbol_t.items():
        for c in j:
            if c in data.keys():
                un_hpo_t[i] = list(set(un_hpo_t[i] + data[c]))
    print(len(un_hpo_t))
    num = 0
    zero_num = []
    for j, i in un_hpo_t.items():
        if len(i) > 0:
            num = num + len(i)
        else:
            zero_num.append(j)
    print(num)
    print(zero_num)
    print(len(zero_num))
    with open('../data/disgenet/DisGeNet_uniprot_hpo.json', 'w') as fggg:
        json.dump(un_hpo_t, fggg, indent=2)


def new():
    # data = pd.read_csv('gene_associations.tsv', sep='\t')
    # print(data[['geneSymbol']].values)
    # gene_symbol = []
    # for i in data[['geneSymbol']].values:
    #     gene_symbol.append(i[0])
    # print(len(list(set(gene_symbol))))
    symbol_hpo = defaultdict(list)
    no_list = []
    with open('gene_symbol_21666.json', 'r') as f:
        data = json.load(f)
    with open('DisGeNet_disease_hpo.json', 'r') as fg:
        DisGeNet_disease_hpo = json.load(fg)
    for i in tqdm(data):
        path = 'gene_disease_total/{}.csv'.format(i)
        if os.path.exists(path):
            data1 = pd.read_csv(path, sep='\t')
        hpo_list = []
        # print(data1.values)
        for j in data1.values:
            # print(j)
            key = j[0].replace('"', '').split(',')[1]
            # print(key)
            if key in DisGeNet_disease_hpo.keys():
                hpo_list = list(set(hpo_list + DisGeNet_disease_hpo[key]))

        if len(hpo_list) > 0:
            symbol_hpo[i] = symbol_hpo[i] + hpo_list
        else:
            no_list.append(i)
    # with open('no_list_total.json', 'w') as fgg:
    #     json.dump(no_list, fgg, indent=2)
    print(len(symbol_hpo))
    num = 0
    zero_num = []
    has_hpo_gene = []
    use_hpo = []
    for j, i in symbol_hpo.items():
        num = num + len(i)
        use_hpo = list(set(use_hpo + i))
        has_hpo_gene.append(j)
    print(len(use_hpo))
    print(len(has_hpo_gene))
    print(num)
    # with open('DisGeNet_gene_hpo_total.json', 'w') as fggg:
    #     json.dump(symbol_hpo, fggg, indent=2)
    # with open('has_hpo_gene.json', 'w') as fgggg:
    #     json.dump(has_hpo_gene, fgggg, indent=2)
    # with open('use_hpo.json', 'w') as fggggg:
    #     json.dump(use_hpo, fggggg, indent=2)


def unpriot_hpo():
    data = pd.read_csv('gene_associations.tsv', sep='\t')
    print(data[['geneId', 'geneSymbol']].values)
    gene_uniprot = pd.read_csv('mapa_geneid_4_uniprot_crossref.tsv', sep='\t')
    print(gene_uniprot.values)
    with open('DisGeNet_gene_hpo_total.json', 'r') as fggg:
        DisGeNet_gene_hpo_total = json.load(fggg)
    uniprot_hpo = {}
    gene_id = {}
    uniprot_list = set()
    id_uniprot = defaultdict(list)
    hpo = []
    for i in data[['geneId', 'geneSymbol']].values:
        gene_id[i[1]] = i[0]
    for i in gene_uniprot.values:
        id_uniprot[i[1]].append(i[0])
    print(gene_id)
    num = 0
    for j, e in DisGeNet_gene_hpo_total.items():
        index_id = gene_id[j]
        if index_id in id_uniprot.keys():
            for k in id_uniprot[index_id]:
                uniprot_hpo[k] = e
                uniprot_list.add(k)
                num = num + len(e)
    print(len(uniprot_hpo))
    print(len(list(uniprot_list)))
    print(num)
    for i in uniprot_hpo.values():
        hpo = list(set(hpo + i))
    print(len(hpo))
    # with open('uniprot_hpo_total.json', 'w') as fggg:
    #     json.dump(uniprot_hpo, fggg, indent=2)
    with open('use_hpo_uniprot.json', 'w') as fgggg:
        json.dump(hpo, fgggg, indent=2)
    # with open('use_uniprot.json', 'w') as fggggg:
    #     json.dump(list(uniprot_list), fggggg, indent=2)


if __name__ == '__main__':
    with open('uniprot_hpo_total.json', 'r') as fggg:
        uniprot_hpo_total = json.load(fggg)
    hpo_uniprot = defaultdict(set)
    hpo_1_10 = []
    hpo_11_30 = []
    hpo_31_100 = []
    hpo_101_300 = []
    hpo_301 = []
    hpo = []
    for i, j in uniprot_hpo_total.items():
        hpo = list(set(hpo + j))
        for k in j:
            hpo_uniprot[k].add(i)
    print(len(hpo))
    print(len(hpo_uniprot))
    for a, b in hpo_uniprot.items():
        if 0 < len(list(b)) <= 10:
            hpo_1_10.append(a)
        elif 10 < len(list(b)) <= 30:
            hpo_11_30.append(a)
        elif 30 < len(list(b)) <= 100:
            hpo_31_100.append(a)
        elif 100 < len(list(b)) <= 300:
            hpo_101_300.append(a)
        elif 300 < len(list(b)):
            hpo_301.append(a)
    print(len(hpo_1_10))
    print(len(hpo_11_30))
    print(len(hpo_31_100))
    print(len(hpo_101_300))
    print(len(hpo_301))
    print(len(hpo_301) + len(hpo_101_300) + len(hpo_31_100) + len(hpo_11_30) + len(hpo_1_10))
    with open('hpo_1_10.json', 'w') as fg:
        json.dump(hpo_1_10, fg, indent=2)
    with open('hpo_11_30.json', 'w') as fgg:
        json.dump(hpo_11_30, fgg, indent=2)
    with open('hpo_31_100.json', 'w') as fggg:
        json.dump(hpo_31_100, fggg, indent=2)
    with open('hpo_101_300.json', 'w') as fgggg:
        json.dump(hpo_101_300, fgggg, indent=2)
    with open('hpo_301.json', 'w') as fggggg:
        json.dump(hpo_301, fggggg, indent=2)
    # unpriot_hpo()