#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author = 'liuloolew@163.com'
@时间：202004141100  v0003
@用途：待定
@说明
@更新：
"""
import argparse, sys, os, re,collections, time, glob
from functools import lru_cache

col_name_single = ["条形码", "单位", "姓名", "样品类型", "检测项目", "文库编号", "primer mix", "index",
                   "文库长度（bp）", "文库浓度（ng/ul)", "稀释倍数", "稀释后浓度（ng/ul)", "体积（ul）",
                   "测序数据量（M）", "摩尔浓度（nM）", "pooling体积（ul）", "质量总量（ng）",
                   "摩尔总量（nM）", "对照样品编号"]

col_name_double = ["条形码", "单位", "姓名", "样品类型", "检测项目", "文库编号", "primer mix", "index",
                   "i5 index", "i7 index", "文库长度（bp）", "文库浓度（ng/ul)", "稀释倍数",
                   "稀释后浓度（ng/ul)", "体积（ul）", "测序数据量（M）", "摩尔浓度（nM）",
                   "pooling体积（ul）", "质量总量（ng）", "摩尔总量（nM）", "对照样品编号"]

primer_mix_indexes_file = "/data_center_17/Project/v1.2-illumina/bin/primer_mix_index_seq.txt"

class Params(object):
    @property
    @lru_cache(1)
    def parser(self):
        return argparse.ArgumentParser(prog='modify_pooling.py',
                                       description=__doc__,
                                       formatter_class=argparse.RawDescriptionHelpFormatter,
                                       epilog="有需要请联系作者：liuloolew@163.com")

    def add_argument(self):
        self.parser.add_argument('-i', '--input_file', dest='input_file',
                                 metavar='FILE', type=str, required=True,
                                 help="pooling文件")
        self.parser.add_argument('-6', '--six_adapter', dest='six_adapter',
                                 metavar='SIXBASES', type=str, default='GATCGGAAGAGCACACGTCT',
                                 help="index长度为6的接头的序列")
        self.parser.add_argument('-8', '--eight_adapter', dest='eight_adapter',
                                 metavar='EIGHTBASES', type=str, default='CTGTCTCTTATACACATCTC',
                                 help="index长度为8的接头的序列")
        self.parser.add_argument('-10', '--ten_adapter', dest='ten_adapter',
                                 metavar='TENBASES', type=str, default='AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA',
                                 help="index长度为10的接头的序列，华大测序仪")
        self.parser.add_argument('-o', '--outdir', dest='outdir',
                                 metavar='DIR', type=str, default='./',
                                 help="输出结果文件夹的名称")
        self.parser.add_argument('-t', '--experiment_time', dest='experiment_time',
                                 metavar='INT', type=str, default=time.strftime("%Y%m%d"),
                                 help="设置实验时间，这里默认设置为当前时间，真正的项目需要写成上机时间，如：20200413")
        self.parser.add_argument('--check', dest='check_first_line',
                                 action='store_true',
                                 help="检查配置文件的行首[列名]与脚本预测的是否一致，主要检测规范的商业样本的命名，对测试样本会报错")
        self.parser.add_argument('-S', '--single_index', dest='single_index',
                                 action='store_true', default=False,
                                 help="是否是单个index测序，模式是单端")
        self.parser.add_argument('-r', '--reverse_seq',dest='reverse_seq', action='store_true',default=False,
                                 help="是否反向互补index序列，默认是False")
        self.parser.add_argument('--test', dest='test_if_or_not', action='store_true',
                                 help="是否为测试样本，若是测试样本则，不会检查样本的index，对照，等等只会按照正常输出")
        self.parser.add_argument('--primer_mix_indexes', dest='primer_mix_indexes',
                                 metavar='FILE', type=str, default=primer_mix_indexes_file,
                                 help="primer_mix对用index的序列文件")
        self.parser.add_argument('--raw_dir', dest='raw_dir',
                                 metavar='DIR', type=str,
                                 help="武汉、南京等通过云上传输的下机数据下载路径，必须是已经合并后的文件")
    def read_params(self):
        self.add_argument()
        args = self.parser.parse_args()
        return args

    @property
    def col_name_single(self):
        if self.params.single_index:
            return col_name_single
        else:
            return col_name_double
        # return col_name_single if self.params.single_index else col_name_double

    # @property
    # def col_name_double(self):
    #     return col_name_double

    def check_primer_mix_indexes(self, sample_name, i5i7_seq_tuple, primer_mix_indexes_file):
        if len(i5i7_seq_tuple) == 3:
            (primer_mix, i5_index, i7_index) = i5i7_seq_tuple
        else:
            (primer_mix, i7_index) = i5i7_seq_tuple

        if primer_mix not in primer_mix_indexes_file:
            print("\t注意: %s primer mix 为 %s, 不在序列信息表格中，请联系项目管理。" % (sample_name, primer_mix))
            return False
        else:
            if len(i5i7_seq_tuple) == 3:
                tmm = i5_index+"_"+i7_index
            else:
                tmm = "/_"+i7_index

            if (tmm != primer_mix_indexes_file[primer_mix]):
                print("\t注意: %s index序列【%s】与序列信息表格【%s】不一致，请联系项目管理。"
                      % (sample_name, tmm, primer_mix_indexes_file[primer_mix]))
                return False
            else:
                return True

    def read_primer_mix_indexes_file(self, primer_mix_indexes):
        dic_ = dict()
        with open(primer_mix_indexes) as inf:
            first_line = inf.__next__()
            for lines in inf:
                tabs = [i.strip() for i in lines.strip().split('\t')]
                index = "_".join(tabs[1:])
                dic_[tabs[0]] = index
        return dic_

    def judge_first_line(self, first_line, col_name):
        first_line_list = first_line.strip().split("\t")
        first_line_list = [i.strip() for i in first_line_list]
        for i, v in enumerate(first_line_list):
            if v != col_name[i]:
                return False
        return True

    def check_line(self, tabs_list):
        flag = "A-I"
        if 'RNA' in tabs_list[1].upper():
            flag = 'R-' if self.params.single_index else 'R-'
        elif 'DNA' in tabs_list[1].upper():
            flag = 'D-' if self.params.single_index else 'D-'

        id_ = flag.join((tabs_list[0], str(tabs_list[3])))
        if "U" in str(tabs_list[3]):
            if self.params.reverse_seq:
                print("\t注意: %s primer mix 为 U 开头, 运行脚本不需要加 -r ,不进行index反向互补,请知晓。" % id_)
        elif "I" in str(tabs_list[3]):
            if self.params.reverse_seq:
                print("\t注意: %s primer mix 为 I 开头, 运行脚本不需要加 -r ,不进行index反向互补,请知晓。" % id_)
        elif "B" in str(tabs_list[3]):
            if not self.params.reverse_seq:
                print("\t注意: %s primer mix 为 B 开头, 运行脚本需要应该加 -r 进行i5反向互补，i7不变,请知晓。" % id_)
        elif "P" in str(tabs_list[3]):
            if not self.params.reverse_seq:
                print("\t注意: %s primer mix 为 P 开头, 运行脚本需要应该加 -r 进行i5和i7反向互补,请知晓。" % id_)
        elif "M" or "A" in str(tabs_list[3]):
            # 武汉,南京数据不需要判断
            pass
        else:
            print("\t注意: %s primer mix 不为 U|I|B|P 开头, 请询问项目管理。" % id_)
        if tabs_list[2] == id_:
            return False
        else:
            return True

    def adapt_info(self, index, six_adpter, night_adpter, ten_adapter):
        if len(index) == 6:
            return six_adpter
        elif len(index) == 8:
            return night_adpter
        elif len(index) == 10:
            return ten_adapter
        elif ten_adapter:
            return ten_adapter

    def which_index_(self, tabs, adapt_):
        return [tabs[4], adapt_, tabs[5]] if self.params.single_index else [tabs[4], tabs[5], adapt_, tabs[6]]

    def get_rev_seq(self, seq):
        return  seq[::-1].translate(str.maketrans('ACGTacgt', 'TGCAtgca'))

    def get_pooling_info(self, params, col_name):
        dic_info = dict()
        flag = 0
        controlsamplename = list()
        sampleinfo_col_ = re.sub(',|，', ',', params.sampleinfo_col)
        index_col = sampleinfo_col_.strip().split(",")
        index_col = [int(i.strip()) for i in index_col]
        # i5和i7是否重复,判断i5和i7对应的prime mix名称是否一致
        i5_i7_index_seq = list()
        primer_mix_indexes_dic = self.read_primer_mix_indexes_file(params.primer_mix_indexes)

        with open(params.input_file) as inf:
            first_line = inf.__next__()
            # 根据参数，是否检查首行
            if params.check_first_line:
                if not self.judge_first_line(first_line, col_name):
                    print("文件列名称：", first_line)
                    print("预测列名称：", col_name)
                    raise Exception('\t 列名称与脚本预设不合适\n')

            for lines in inf:
                if not len(lines.strip()) or lines.strip().startswith('#'):  # 判断是否是空行或注释行
                    continue
                tabs = [i.strip() for i in lines.strip().split('\t')]
                tabs = [tabs[i] for i in index_col]
                if self.check_line(tabs, ):
                    raise Exception('\t\t %s 样本的index编号有问题\n' % tabs[0])


                # 判断index和i7片段是否有重复  双index所在的tab[4:6]  ; 单index tab[4]
                i5i7_seq_tuple = tuple(tabs[3:5]) if params.single_index else tuple(tabs[3:6])
                # 20200828 检擦index名称对用的序列名称 [tabs[3], tabs[4],tabs[5]] U1 GCGCATAT CTGATCGT
                self.check_primer_mix_indexes(tabs[0], i5i7_seq_tuple, primer_mix_indexes_dic)

                if i5i7_seq_tuple in i5_i7_index_seq:
                    print('\t\t %s 样本的i5 index 和 i7 index 有重复\n' % tabs[0])
                else:
                    i5_i7_index_seq.append(i5i7_seq_tuple)

                sample_id = tabs[2]
                if "RNA" in tabs[1]:
                    sample_id = '%s_RNA' % tabs[2]
                elif "DNA" in tabs[1]:
                    sample_id = '%s_DNA' % tabs[2]
                adapt_ = self.adapt_info(tabs[4], params.six_adapter, params.eight_adapter, params.ten_adapter)
                if "NC" in tabs[0]:
                    controlsamplename.append(sample_id)
                elif not tabs[0].startswith('RP'):
                    # raise Exception('\t\t %s 样本名称不是【RP】开头\n' % tabs[0])
                    print('注意:\n\t %s 样本名称不是【RP】开头\n' % tabs[0])

                if "NC" not in tabs[0]:
                    if flag not in dic_info:
                        dic_info[flag] = collections.OrderedDict()
                    if sample_id not in dic_info[flag]:
                        dic_info[flag][sample_id] = self.which_index_(tabs, adapt_)
                    else:
                        # raise Exception('\t\t %s 样本在分析表中有重复\n' % tabs[0])
                        print('\t\t %s 样本在分析表中有重复\n' % tabs[0])
                else:
                    flag += 1
                    if flag not in dic_info:
                        dic_info[flag] = collections.OrderedDict()
                    if sample_id not in dic_info[flag]:
                        dic_info[flag][sample_id] = self.which_index_(tabs, adapt_)  # RP20JSUBR0541D-I16, CCGTCC, 20D21NC3
                    else:
                        # raise Exception('\t\t %s 样本在分析表中有重复\n' % tabs[0])
                        print('\t\t %s 样本在分析表中有重复\n' % tabs[0])
        return dic_info

    def sample_config_file(self, dic_, outdir):
        dic_sample_config_file = collections.OrderedDict()
        # 取出 对照样本的id列表
        controlsamplename_name_list = list(dic_.keys())[1:]
        controlsamplename_name = [list(dic_[i].keys())[0] for i in controlsamplename_name_list]
        for flag in dic_.keys():
            for sample in dic_[flag]:
                # 20200423 增加一个判断，提供的第18列[实验室提供的对照品的名称] 是否和我们对照品一致
                lab_compare = dic_[flag][sample][-1]
                if sample == controlsamplename_name[flag - 1]:
                    script_lab_compare = controlsamplename_name[flag - 1]
                else:
                    script_lab_compare = controlsamplename_name[flag]
                if lab_compare not in script_lab_compare:
                    raise Exception('\t\t %s脚本得到的对照品名称[%s]与实验室提供[%s]的不一致，请检查...\n'
                                    % (sample, script_lab_compare, lab_compare))
                # adapter 碱基名称
                index_base = dic_[flag][sample][-2]
                # 对照品
                if sample == controlsamplename_name[flag - 1]:
                    values_dic_sample = '%s\t%s\tN\t%s' % (sample, controlsamplename_name[flag - 1], index_base)
                    dic_sample_config_file[sample] = values_dic_sample
                else:
                    # 商业样本
                    values_dic_sample = '%s\t%s\tN\t%s' % (sample, controlsamplename_name[flag], index_base)
                    dic_sample_config_file[sample] = values_dic_sample
        with open("%s/sample_config.txt" % outdir, 'w') as out:
            out.write('#sampleid\tcontrolsamplename\tmapunique\tadapter\n')
            for k in dic_sample_config_file:
                out.write('%s\n' % dic_sample_config_file[k])
            return dic_sample_config_file


    def wh_sample_config_file(self, dic_, outdir, raw_dir):
        samples_index_num_dic = dict()
        for file in glob.glob('%s/*.fq.gz' % raw_dir):
            try:
                m = re.search(r'_(\d+).fq.gz', file).group(1)  # /data_center_16/DNA_Data/realpatho_wuhan/BN20H20S75-1/V300071801_94.fq.gz
                samples_index_num_dic[m] = file
            except:
                # /data_center_16/DNA_Data/realpatho_nanjing/20201010_S200014172_Result/Rawdata/S200014172_L01_77-80.fq.gz
                m = re.search(r'_(\d+-\d+).fq.gz', file).group(1)
                samples_index_num_dic[m] = file

        with open("%s/sample_config.txt" % outdir, "a") as wh_out:
            wh_out.write("#add\n")
            for sample in dic_:
                try:
                    # 210506 修改正则匹配
                    try:
                        m = re.search(r'-([a-zA-Z]?)(\w+-\w+)_', sample).group(2)  # RP20BHUZY0015D-M41_DNA RP20JSUBR2931D-A13-14_DNA
                    except:
                        m = re.search(r'-([a-zA-Z]?)(\d+)_', sample).group(2)  # RP20BHUZY0015D-41_DNA  RP20BHUZY0015D-M41_DNA
                    if m in samples_index_num_dic:
                        wh_out.write("%s\trawdata\t%s\n" % (sample ,samples_index_num_dic[m]))
                    else:
                        print('注意:\n\t %s 样本的index序列名称在rawdata中未找到\n' % sample)
                except:
                    pass


    def detailed_information(self, experiment_time):
        year, mouth, day = experiment_time[:4], experiment_time[4:6], experiment_time[6:]
        detailed_information_list = list()
        detailed_information_list.append('[Header],,,,,')
        detailed_information_list.append('Local Run Manager Analysis Id,101,,,,')
        detailed_information_list.append('Date,%d/%d/%s,,,,' % (int(mouth), int(day), year))
        detailed_information_list.append('Experiment Name,%s,,,,' % experiment_time)
        detailed_information_list.append('Workflow,GenerateFastQWorkflow,,,,')
        detailed_information_list.append(
            'Description,Auto generated sample sheet.  Used by workflow module to kick off Isis analysis,,,,')
        detailed_information_list.append('Chemistry,Amplicon,,,,\n,,,,,')
        detailed_information_list.append('[Reads],,,,,\n75,,,,,\n75,,,,,\n,,,,,')
        detailed_information_list.append('[Settings],,,,,\nAdapter,,,,,\n,,,,,\n[Data],,,,,')
        detailed_information_list.append('Sample_ID,Sample_Name,index,I7_Index_ID,index2,I5_Index_ID')
        return detailed_information_list

    def get_caifen(self, dic_, outdir, experiment_time):
        detailed_information_list = self.detailed_information(experiment_time)

        split_dict = dict()
        for flag in dic_.keys():
            for sample in dic_[flag]:
                index_num_= 1 if self.params.single_index else 2
                index_base = dic_[flag][sample][:index_num_]
                index_base_len = len(index_base[0])

                if index_base_len not in split_dict:
                    split_dict[index_base_len] = dict()
                split_dict[index_base_len][sample] = index_base

        for index_base_len in split_dict:
            i = list(split_dict.keys()).index(index_base_len)
            with open("%s/SampleSheetUsed.%d.csv" % (outdir, i + 1), 'w') as  out2:
                out2.write("%s\n" % '\n'.join(detailed_information_list))
                for sample in split_dict[index_base_len]:
                    # I开头,单端不需要反转
                    if self.params.single_index:
                        out2.write('%s,%s,%s,,,\n' % (sample, sample, split_dict[index_base_len][sample][0]))
                    else:
                        # index 反转互补
                        if self.params.reverse_seq:
                            # P端是i7、i5都反向互补；RP20SDUBR0362D-U24_DNA
                            if "-P" in sample[-8:-4]:
                                out2.write('%s,%s,%s,,%s,\n' % (sample, sample,
                                                                self.get_rev_seq(split_dict[index_base_len][sample][1]),
                                                                self.get_rev_seq(split_dict[index_base_len][sample][0])))
                            # B端是i7不反向互补，i5反向互补。
                            elif "-B" in sample[-8:-4]:
                                out2.write('%s,%s,%s,,%s,\n' % (sample, sample, split_dict[index_base_len][sample][1],
                                                                self.get_rev_seq(split_dict[index_base_len][sample][0])))
                        # U端是i7、i5都不反向互补；
                        else:
                            out2.write('%s,%s,%s,,%s,\n' % (sample, sample, split_dict[index_base_len][sample][1],
                                                            split_dict[index_base_len][sample][0]))


    def __call__(self):
        self.params = self.read_params()

        if not self.params.test_if_or_not:
            sample_info_dic = self.get_pooling_info(self.params, self.col_name_single)
            sample_config_file_dic = self.sample_config_file(sample_info_dic, self.params.outdir)
            if self.params.raw_dir:
                self.wh_sample_config_file(sample_config_file_dic, self.params.outdir, self.params.raw_dir)
            else:
                self.get_caifen(sample_info_dic, self.params.outdir, self.params.experiment_time)


class SingleIndex(Params):
    def __init__(self):
        super().__init__()
#
    def add_argument(self):
        super(SingleIndex, self).add_argument()
        self.parser.add_argument('-s', '--sampleinfo_col', dest='sampleinfo_col',
                                 metavar='STR', type=str, default='0,4,5,6,7,18',
                                 help="构建参数文件所需的列的编号(根据pooling文件)，分别是"
                                      "[条形码, 检测项目, 文库编号, index, index序列, 对照样品编号]")


if __name__ == '__main__':
    run = SingleIndex()
    run()
