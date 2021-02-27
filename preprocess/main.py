import numpy
import os


# 处理原始临床数据 将样本id存入字典
def processing_clinical_data( ):

    pdata_in_clinical_dict_all = {}
    pdata_in_clinical_dict = {}
    label_dict = {}
    n = 0
    with open(raw_clinical_file) as f:
        for l in f.readlines():
            if len(l) > 0:
                l = l.strip('\n').split('\t')
                data = [float(i) for i in l[1:]]
                name = l[0]
                n = n + 1
                pdata_in_clinical_dict_all[name] = data

    clinical_patients_num = 0

    for key in pdata_in_clinical_dict_all:
        data = pdata_in_clinical_dict_all[key]
        clinical_patients_num += 1
        n = 1
        dd = []
        for i in range(len(data)):
            if n == label_position:
                if data[i] > 0:
                    label_dict[key] = data[i]
            else:
                dd.append(data[i])
            n += 1

        pdata_in_clinical_dict[key] = dd

    print('label_dict:')
    print(label_dict)
    print('\n')

    print('pdata_in_cinical_dict:')
    print(pdata_in_clinical_dict)
    print('\n')

    print('patients number in clinical: ', clinical_patients_num)

    pdata_in_cinical_file = pathdict + '/pdata_in_clinical_dict.txt'
    cid = open(pdata_in_cinical_file, 'w')
    cid.write(str(pdata_in_clinical_dict))

    label_file = pathdict + '/label_dict.txt'
    cid = open(label_file, 'w')
    cid.write(str(label_dict))


# 处理原始基因表达数据 将原始基因表达数据转换 每一行对应一个样本 与临床数据一一对应
def processing_gene_cna():

    gene_file = pathother + '/data_gene_unmatch.txt'
    cna_file = pathother + '/data_cna_unmatch.txt'

    pid_in_gen_dict = {}
    gid_in_gene_dict = {}

    pid_in_cna_dict = {}
    gid_in_cna_dict = {}

    gene_wf = open(gene_file, 'w')
    cna_wf = open(cna_file, 'w')

    n = 0
    with open(raw_gene_file) as f:
        for l in f.readlines():
            if len(l) > 0:
                d = l.strip('\n').split('\t')
                pid = d[0]
                data = [str(i) for i in d[1:]]
                if n > 0:
                    gid_in_gene_dict[pid] = n
                    for i in range(len(data)):
                        if data[i] == 'NA':
                            data[i] = 0
                        # print(n,i)
                        if (i + 1 == len(data)):
                            if float(data[i]) < gene_min:
                                gene_wf.write('-1\n')
                            elif float(data[i]) > gene_max:
                                gene_wf.write('1\n')
                            else:
                                gene_wf.write('0\n')
                        else:
                            if float(data[i]) < gene_min:
                                gene_wf.write('-1 ')
                            elif float(data[i]) > gene_max:
                                gene_wf.write('1 ')
                            else:
                                gene_wf.write('0 ')
                if n == 0:
                    for i in range(len(data)):
                        pid_in_gen_dict[data[i]] = i
                n = n + 1
    k = 0
    with open(raw_cna_file) as cnaf:
        for line in cnaf.readlines():
            if len(line) > 0:
                d = line.strip('\n').split('\t')
                cna = d[0]
                data = [str(i) for i in d[1:]]

                if k > 0:
                    gid_in_cna_dict[cna] = k
                    for i in range(len(data)):
                        if data[i] == 'NA':
                            data[i] = 0
                        if i + 1 == len(data):
                            cna_wf.write(str(data[i]) + '\n')
                        else:
                            cna_wf.write(str(data[i]) + ' ')
                if k == 0:
                    for i in range(len(data)):
                        pid_in_cna_dict[data[i]] = i
                k = k + 1
    gid_gen = pathdict + '/gid_in_gene_dict.txt'
    pid_gen = pathdict + '/pid_in_gene_dict.txt'
    gid_cna = pathdict + '/gid_in_cna_dict.txt'
    pid_cna = pathdict + '/pid_in_cna_dict.txt'

    gid_gen_file = open(gid_gen, 'w')
    pid_gen_file = open(pid_gen, 'w')
    gid_cna_file = open(gid_cna, 'w')
    pid_cna_file = open(pid_cna, 'w')

    gid_gen_file.write(str(gid_in_gene_dict))
    pid_gen_file.write(str(pid_in_gen_dict))
    gid_cna_file.write(str(gid_in_cna_dict))
    pid_cna_file.write(str(pid_in_cna_dict))


def matrix():
    cna_matrix_file = pathother + '/data_cna_unmatch.txt'
    cna_matrix = numpy.loadtxt(cna_matrix_file, delimiter=' ')
    c_matrix = numpy.transpose(cna_matrix)

    gene_matrix_file = pathother + '/data_gene_unmatch.txt'
    gene_matrix = numpy.loadtxt(gene_matrix_file, delimiter=' ')
    g_matrix = numpy.transpose(gene_matrix)

    # print('gene matrix:', g_matrix.shape)
    gene_num = g_matrix.shape[1]
    print('patients number in gene:', g_matrix.shape[0])
    print('gene number:', g_matrix.shape[1])
    print('\n')
    # print('cna matrix:', c_matrix.shape)
    cna_num = c_matrix.shape[1]
    print('patients number in cna:', c_matrix.shape[0])
    print('cna number:', cna_num)

    # numpy.savetxt(pathdata + '/data_cna_match_final.txt', c_matrix, fmt='%.0f')
    # numpy.savetxt(pathdata + '/data_gene_match_final.txt', g_matrix, fmt='%.0f')

    gene_file = pathdata + '/%s_gene.txt' % (data_des)
    cna_file = pathdata + '/%s_cna.txt' %(data_des)
    clinical_file = pathdata + '/%s_clinical.txt' % (data_des)
    label_file = pathother + '/label_real_number.txt'

    gene_wf = open(gene_file, 'w')
    cna_wf = open(cna_file, 'w')
    clinical_wf = open(clinical_file, 'w')
    label_wf = open(label_file, 'w')

    f1 = open(pathdict + '/pid_in_gene_dict.txt', 'r')
    f2 = open(pathdict + '/pid_in_cna_dict.txt', 'r')
    f3 = open(pathdict + '/pdata_in_clinical_dict.txt', 'r')
    f4 = open(pathdict + '/label_dict.txt', 'r')

    pid_in_gene_dict = eval(f1.read())
    pid_in_cna_dict = eval(f2.read())
    pdata_in_clinical_dict = eval(f3.read())
    label_dict = eval(f4.read())

    pdata_in_gene_dict = {}
    pdata_in_cna_dict = {}

    i = 0
    for key in pid_in_gene_dict:
        pdata_in_gene_dict[key] = g_matrix[i]
        i = i + 1
    k = 0
    for key in pid_in_cna_dict:
        pdata_in_cna_dict[key] = c_matrix[k]
        k = k + 1
    patient_number_final = 0
    for key in pdata_in_clinical_dict:
        if key in label_dict:
            if key in pdata_in_gene_dict:
                if key in pdata_in_cna_dict:

                    patient_number_final += 1

                    label_wf.write(str(label_dict[key]))
                    label_wf.write('\n')

                    clinical = pdata_in_clinical_dict[key]
                    for i in range(len(clinical)):
                        if i + 1 == len(clinical):
                            clinical_wf.write(str(clinical[i]) + '\n')
                        else:
                            clinical_wf.write(str(clinical[i]) + ' ')

                    gene = pdata_in_gene_dict[key]
                    cna = pdata_in_cna_dict[key]
                    for i in range(len(gene)):
                        if i + 1 == len(gene):
                            gene_wf.write(str(gene[i]) + '\n')
                        else:
                            gene_wf.write(str(gene[i]) + ' ')
                    for i in range(len(cna)):
                        if i + 1 == len(cna):
                            cna_wf.write(str(cna[i]) + '\n')
                        else:
                            cna_wf.write(str(cna[i]) + ' ')
    f1.close()
    f2.close()
    f3.close()
    gene_wf.close()
    cna_wf.close()

    new_path = './data/%s_pnum%s_gnum%s_cnum%s' % (data_des, patient_number_final, gene_num, cna_num)
    try:
        os.rename(pathdata, new_path)
    except Exception as e:
        print(e)
        print('rename dir fail\r\n')
    else:
        print('rename dir success\r\n')
    return new_path

def label():
    label_file = new_path + '/other/label_real_number.txt'
    label_final_file = new_path + '/%s_label.txt'%(data_des)
    label_wf = open(label_final_file, 'w')

    pos = 0
    neg = 0
    with open(label_file) as fl:
        for ll in fl.readlines():
            data = float(ll)
            if data > label_cut:
                label_wf.write('1\n')
                pos += 1
            else:
                label_wf.write('0\n')
                neg += 1
    print("pos:", pos)
    print("neg:", neg)

    new_label_path  = new_path + '/%s_label_cut%sm_p%s_n%s.txt'%(data_des, label_cut, pos, neg)
    try:
        os.rename(label_final_file, new_label_path)
    except Exception as e:
        print(e)
        print('rename dir fail\r\n')
    else:
        print('rename dir success\r\n')

def train_datatype():
    gen_file = new_path + '/%s_gene.txt'%(data_des)
    cna_file = new_path + '/%s_cna.txt'%(data_des)
    train_gene_file = new_path + '/%s_train_gene.txt'%(data_des)
    train_cna_file = new_path + '/%s_train_cna.txt'%(data_des)
    rf1 = open(train_gene_file, 'w')
    rf2 = open(train_cna_file, 'w')

    k = 0
    with open(gen_file) as f:
        for l in f.readlines():
            l = l.strip('\n').split(' ')
            if len(l) > 0:
                rf1.write(str(k))
                for i in range(len(l)):
                    if float(l[i]) > 0.0 or float(l[i]) < 0.0:
                        rf1.write(' ' + str(i))
                rf1.write("\n")
            k += 1
    n = 0
    with open(cna_file) as ff:
        for l in ff.readlines():
            l = l.strip('\n').split(' ')
            if len(l) > 0:
                rf2.write(str(n))
                for i in range(len(l)):
                    if float(l[i]) > 0.0 or float(l[i]) < 0.0:
                        rf2.write(' ' + str(i))
                rf2.write("\n")
            n += 1


data_path = ''
gene_max = 1.0
gene_min = -1.0
label_cut = 12 * 2
label_position = 25
#colorectal 22
#breast
#lung
#lymphoid 32   all 41
#breast actg 27
data_des = 'breast_cancer_data_acgt_clin1'

raw_clinical_file = './raw_data_breast_tcga_clin/clinical.txt'
raw_cna_file = './raw_data_breast_tcga_clin/cna.txt'
raw_gene_file = './raw_data_breast_tcga_clin/gene.txt'

pathdata = './data/%s' % (data_des)
if not os.path.exists(pathdata):
    os.mkdir(pathdata)

pathdict = './data/%s/dict' % (data_des)
if not os.path.exists(pathdict):
    os.mkdir(pathdict)

pathother = './data/%s/other' % (data_des)
if not os.path.exists(pathother):
    os.mkdir(pathother)


processing_clinical_data()

processing_gene_cna()

new_path = matrix()

label()

train_datatype()

# train_or_test(pathdata)
