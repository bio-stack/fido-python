from lxml.etree import iterparse
import string
import sys
import copy
import networkx as nx

import pylab as P

import psyco
psyco.full()


NS_prefix='{http://regis-web.systemsbiology.net/pepXML}'
PEPTIDE_THRESHOLD = 0.2

def cleanupPeptide(pepStr):
    return pepStr.replace('I','L')

def proteinPeptideEdges(filename):
    result = set()

    namespaces = []
    for (event, node) in iterparse(filename, ['start', 'end', 'start-ns', 'end-ns']):
        if event == 'start-ns':
            namespaces.insert(0, node)
        elif event == 'end-ns':
            namespaces.pop(0)
        elif event == 'end':
            namespace, tag = string.split(node.tag[1:], '}', 1)
            if tag == 'search_hit':
                # make sure this peptide passes the score threshold
                peptideprophet_result = node.find(NS_prefix + "analysis_result/" + NS_prefix + "peptideprophet_result")

                if peptideprophet_result != None:
                    peptideProbability = float(peptideprophet_result.attrib.get('probability'))
                    if peptideProbability < PEPTIDE_THRESHOLD:
                        continue
                elif peptideProbability < PEPTIDE_THRESHOLD:
                    continue

                # get the protein and peptide from this match
                protein = node.attrib.get("protein")
                peptide = node.attrib.get("peptide")

                result.add((protein, peptide))

                alternative_proteins = node.findall(NS_prefix + "alternative_protein")
                for alt_protein_node in alternative_proteins:
#                    alt_namespace, alt_tag = string.split(alt_protein_node.tag[1:], '}', 1)
                    alt_protein = alt_protein_node.attrib.get("protein")
#                    print "alt",alt_protein
                    result.add((alt_protein, peptide))

                # process the alternative proteins
#                 for child in node.getchildren():
#                     child_namespace, child_tag = string.split(child.tag[1:], '}', 1)
#                     if child_tag == 'alternative_protein':
#                         alt_protein = child.attrib.get("protein")
#                         result.add((alt_protein, peptide))


    return result

def proteinPeptideEdgesFromPivdo(filename):
    result = set()

    file = open(filename, 'r')

    state = 'e'

    lineNumber = 0
    for line in file.readlines():
        lineNumber += 1
#            print lineNumber
        code,data = line.split()
            
        # some compatability with .pivdo2
        if code == 'd' or code =='c':
            continue

        # code is the prefix for the line
        # state is the next code symbol that the FSM wants
        if code not in state:
            if code == 'e' and 'p' in state:
                print >> sys.stderr, 'Warning: peptide',pepString,'had no probability, using the weight from the previous peptide'
#                    self.addMax(self.peptideDict, pepString, pepVal)
            else:
                print >> sys.stderr, 'Error: on line',lineNumber,'FSM received code',code,'but was expecting code',state
                exit(1)

        if code == 'e':
            pepString = data
            state = 'r'
        elif code == 'r':
            protString = data
            result.add((protString, pepString))
            state = ('r','p')
        elif code in 'p':
            pepVal = float(data)
            state = 'e'
        else:
            print >> sys.stderr, 'Error: on line',lineNumber,'unrecognized FSM code',code

    return result



#res = proteinPeptideEdges('catted_search_results.pepXML')

def readSubgraphsFromPivdo(filename):
    res = proteinPeptideEdgesFromPivdo(filename)
    print "Read", len(res), "entries"
    G = nx.Graph()
    for prot,pep in res:
        pep = cleanupPeptide(pep)
        G.add_edge("R_"+prot, "E_"+pep)
        
    subgraphs = nx.connected_component_subgraphs(G)
    return subgraphs


def count_prots(subgraph):
	prot_count = 0
	for n in subgraph.nodes():
		if n.startswith("R_"):
			prot_count += 1

	return prot_count

def getBigSubgraphs(subgraphs):
    BIG = 10
    bigSG = []
    for sg in subgraphs:
        if count_prots(sg) > BIG:
            bigSG.append(sg)
            #        drawGraph(bigSG)
            #    	nx.draw_graphviz(sg)
            #        P.show()
    return bigSG

def drawGraph(sg):
    protSet = set()
    for ind, node in enumerate(sg.nodes()):
	if node.startswith("R_"):
            protSet.add(ind)
    
    sizes = []
    colors = []
    for i in xrange(0,len(sg)):
	if i in protSet:
            sizes.append(200)
            colors.append('r')
	else:
            sizes.append(50)
            colors.append('b')

    nx.draw_graphviz(sg, node_size = sizes, with_labels = False, node_color = colors)

def getProteinGroupGraph(sg):
    edges = {}
    for n in sg:
	edges[ frozenset([ i for i in sg[n] ]) ] = n
    
    newSG = copy.deepcopy(sg)
    for i in set([i for i in sg]) - set([edges[e] for e in edges]):
	newSG.remove_node(i)

    return newSG

def getProteinToProteinGraph(sg):
    G = nx.Graph()
    for ind, n in enumerate(sg.nodes()):
	if not n.startswith("R_"):
            continue
	for n2 in sg.nodes()[ind+1:]:
            if len(set(sg[n]).intersection(set(sg[n2]))) > 0:
#                print "OVERLAP"
#                print "\t",set(sg[n])
#                print "\t",set(sg[n2])
                G.add_edge(n, n2)
    return G

# for sg in bigSG:
# 	temp = getProteinGroupGraph(sg)
# 	print "ORIGINAL",count_prots(sg)
# 	drawGraph(sg)
# 	P.show()
# 	print "NEW",count_prots(temp)
# 	drawGraph(temp)
# 	P.show()
