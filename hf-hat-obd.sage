## Copyright (C) 2025 Cagatay Kutluhan, Gordana Matic, Jeremy Van Horn-Morris, Andy Wand
## Contact: kutluhan@buffalo.edu
## Based on hf-hat
## Copyright (C) 2015 Sucharit Sarkar
## Source: https://github.com/sucharit/hf-hat
## Contact: sucharit@math.princeton.edu

## hf-hat-obd is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of
## the License, or (at your option) any later version.

## hf-hat-obd is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with hf-hat-obd; see COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

class HeegaardDiagram():

    def __init__(self,boundary_intersections,num_pointed_regions,name=None):
        #Basic initialization. (Will convert boundary_intersections to its lexicographically smallest form.)
        self.name=str(name)
        self.boundary_intersections=[]
        for intersections_list in boundary_intersections:
            smallest_indices=[2*ind_i for ind_i,i in enumerate(intersections_list[0::2]) if i==min(intersections_list[0::2])]
            next_intersections=[intersections_list[ind_i+1] for ind_i in smallest_indices]
            smallest_index=next(ind_i for ind_i in smallest_indices if intersections_list[ind_i+1]==min(next_intersections))
            self.boundary_intersections.append(intersections_list[smallest_index:]+intersections_list[:smallest_index])

        self.regions=range(len(self.boundary_intersections))#the regions
        self.regions_un=range(len(self.boundary_intersections)-num_pointed_regions)#the unpointed regions
        self.intersections=list(range(1+max(flatten(self.boundary_intersections))))#the intersection points

        #Euler measures of the regions, times 2 (unpointed data stored as a matrix).
        self.euler_measures_2=[2-len(self.boundary_intersections[R])/2 for R in self.regions]
        self.euler_measures_2_un=vector(ZZ,len(self.regions_un),self.euler_measures_2[:len(self.regions_un)])
        if min(self.euler_measures_2_un)>-1:
            self.is_nice=True#is a nice diagram
        else:
            self.is_nice=False
        #boundary_mat[i][j] is the coefficient of the boundary of (boundary of the i-th region restricted to alpha circles) at the j-th intersection point  (unpointed data stored as a matrix).
        self.boundary_mat=[[-self.boundary_intersections[R][0::2].count(p)+self.boundary_intersections[R][1::2].count(p) for p in self.intersections] for R in self.regions]
        self.boundary_mat_un=matrix(ZZ,len(self.regions_un),len(self.intersections),self.boundary_mat[:len(self.regions_un)])
        
        #point_measures_4[i][j] is the point measure of i-th region at the j-th point, times 4  (unpointed data stored as a matrix).
        self.point_measures_4=[[self.boundary_intersections[R].count(p) for p in self.intersections] for R in self.regions]
        self.point_measures_4_un=matrix(ZZ,len(self.regions_un),len(self.intersections),self.point_measures_4[:len(self.regions_un)])
        
        self.periodic_domains=kernel(self.boundary_mat_un)#the free Abelian group of periodic domains
        self.pd_basis=basis(self.periodic_domains)#basis in echelon form for periodic domains
        self.b1=len(self.pd_basis)
        
        #intersections_on_alphas[i] is the ordered list of intersections on alpha_i (according to the orientation of alpha_i). Similarly, for beta.
        self.intersections_on_alphas=[]
        while len(flatten(self.intersections_on_alphas))<len(self.intersections):
            start_p=next(p for p in self.intersections if p not in flatten(self.intersections_on_alphas))
            new_circle=[]
            curr_p=start_p
            while curr_p!=start_p or new_circle==[]:
                new_circle.append(curr_p)
                found_next_point=False
                regions_with_curr_point=[(R,ind_p) for R in self.regions for ind_p in range(len(self.boundary_intersections[R])) if curr_p == self.boundary_intersections[R][ind_p]]
                for (R,ind_p) in regions_with_curr_point:
                    if not found_next_point:
                        if ind_p%2==0 and self.boundary_intersections[R][ind_p+1] not in new_circle:
                            found_next_point=True
                            curr_p=self.boundary_intersections[R][ind_p+1]
                        elif ind_p%2==1 and self.boundary_intersections[R][ind_p-1] not in new_circle:
                            found_next_point=True
                            curr_p=self.boundary_intersections[R][ind_p-1]
                if not found_next_point:#must have completed the cycle
                    curr_p=start_p
            self.intersections_on_alphas.append(new_circle)
            
        self.alphas=range(len(self.intersections_on_alphas))#the alpha circles
        self.contact_intersections=[]
        for a in self.alphas:
            self.contact_intersections.append(self.intersections_on_alphas[a][0])
            
        self.intersections_on_betas=[]
        while len(flatten(self.intersections_on_betas))<len(self.intersections):
            start_p=next(p for p in self.contact_intersections if p not in flatten(self.intersections_on_betas))
            new_circle=[]
            curr_p=start_p
            while curr_p!=start_p or new_circle==[]:
                new_circle.append(curr_p)
                found_next_point=False
                regions_with_curr_point=[(R,ind_p) for R in self.regions for ind_p in range(len(self.boundary_intersections[R])) if curr_p == self.boundary_intersections[R][ind_p]]
                for (R,ind_p) in regions_with_curr_point:
                    if not found_next_point:
                        if ind_p%2==0 and self.boundary_intersections[R][ind_p-1] not in new_circle:
                            found_next_point=True
                            curr_p=self.boundary_intersections[R][ind_p-1]
                        elif ind_p%2==1 and self.boundary_intersections[R][ind_p+1-len(self.boundary_intersections[R])] not in new_circle:
                            found_next_point=True
                            curr_p=self.boundary_intersections[R][ind_p+1-len(self.boundary_intersections[R])]
                if not found_next_point:#must have completed the cycle
                    curr_p=start_p
            self.intersections_on_betas.append(new_circle)

        self.betas=range(len(self.intersections_on_betas))#the beta circles
        if sorted(flatten(self.intersections_on_alphas))!=self.intersections:
            raise Exception("Alpha circles don't contain all the intersections")
        if sorted(flatten(self.intersections_on_betas))!=self.intersections:
            raise Exception("Beta circles don't contain all the intersections")
        
        self.adjacent_region_list=[]#list of unpointed regions adjacent to contact intersection points
        for p in self.contact_intersections:
            for region in self.regions_un:
                if p in self.boundary_intersections[region] and region not in self.adjacent_region_list:
                    self.adjacent_region_list.append(region)
        
        N=len(self.contact_intersections)
        self.regions_by_points=[]#unpointed region lists by contact intersection points
        for i in range(len(self.alphas)):
            regions_by_points_i=[]
            for region in self.adjacent_region_list:
                if self.contact_intersections[i] in self.boundary_intersections[region]:
                    regions_by_points_i.append(region)
            self.regions_by_points.append(regions_by_points_i)
    
        self.region_pairs=[]
        for i in range(len(self.alphas)):
            if len(self.regions_by_points[i])==2:
                self.region_pairs.append(self.regions_by_points[i])
        
        self.indicator_list=[]#the list of periodic domains in pd_basis and their multiplicities in the regions from adjacent_region_list
        for pdom in self.pd_basis:
            for region in self.adjacent_region_list:
                self.indicator_list.append(((self.pd_basis).index(pdom),region,pdom[region]))

        m=len(self.adjacent_region_list)
        n=len(self.pd_basis)
        self.indicator_matrix=matrix(ZZ,m,n)#matrix whose columns are multiplicities of periodic domains in adjacent unpointed regions at all contact intersection points
        for i in range(m):
            for j in range(n):
                self.indicator_matrix[i,j]=self.pd_basis[j][self.adjacent_region_list[i]]

        self.intersection_incidence=[[next(a for a in self.alphas if p in self.intersections_on_alphas[a]),next(b for b in self.betas if p in self.intersections_on_betas[b])] for p in self.intersections]#the i-th intersection point lies in alpha_a and beta_b, and has alpha.beta intersection number n, where [a,b,n]=intersection_incidence[i]
        for p in self.intersections:
            try:
                [a,b]=self.intersection_incidence[p]

                if len(self.intersections_on_alphas[a])>2 and len(self.intersections_on_betas[b])>2:#if any of a or b has 2 or fewer intersections, their orientations are not yet well-defined.
                    [a_pind,b_pind]=[(self.intersections_on_alphas[a]).index(p),(self.intersections_on_betas[b]).index(p)]
                    (R,ind_p)=next((R,ind_p) for R in self.regions for ind_p in range(len(self.boundary_intersections[R])/2) if self.boundary_intersections[R][2*ind_p]==p)
                    prev_b=self.boundary_intersections[R][2*ind_p-1]
                    next_a=self.boundary_intersections[R][2*ind_p+1]
                    intersection_number=-1#if a were oriented from p to next_a, and b oriented from prev_b to p, then the intersection.
                    if self.intersections_on_alphas[a][a_pind-1]==next_a:#a is actually oriented oppositely
                        intersection_number=-intersection_number
                    if self.intersections_on_betas[b][b_pind-1]!=prev_b:#b is actually oriented oppositely
                        intersection_number=-intersection_number
                    (self.intersection_incidence[p]).append(intersection_number)
                    
                elif len(self.intersections_on_alphas[a])==1 or len(self.intersections_on_betas[b])==1:
                    (self.intersection_incidence[p]).append(1)
                else:
                    print("WARNING: The current implementation of intersection numbers doesn't work if some alpha or beta circle has exactly 2 intersections.")
                    raise Exception("Couldn't orient the alpha and beta circles")#Comment out exception depending on how we feel.
            except ValueError:#intersection number at p already determined
                pass
        self.intersection_matrix=matrix(ZZ,len(self.alphas),len(self.betas))#the (a,b) entry is the intersection number between alpha circle a and beta circle b
        for a in self.alphas:
            for b in self.betas:
                self.intersection_matrix[a,b]=sum([n for [i,j,n] in self.intersection_incidence if (a==i and b==j)])
        
        #region_graph_alpha is a graph whose vertices are regions, and edges are alpha arcs separating two regions. Edges labeled by which alpha circle, and the alpha arc (as the first (according to the orientation of that alpha circle) intersection point on that alpha arc). Ditto for beta circles
        self.region_graph_alpha=Graph(len(self.regions),multiedges=True,loops=True)
        for a in self.alphas:
            for ind_p,p in enumerate(self.intersections_on_alphas[a]):
                q=self.intersections_on_alphas[a][ind_p+1-len(self.intersections_on_alphas[a])]
                regions_with_pq_1=[R for R in self.regions for ind_foo,foo in enumerate(self.boundary_intersections[R][0::2]) if [foo,self.boundary_intersections[R][2*ind_foo+1]]==[p,q]]
                regions_with_pq_2=[R for R in self.regions for ind_foo,foo in enumerate(self.boundary_intersections[R][0::2]) if [foo,self.boundary_intersections[R][2*ind_foo+1]]==[q,p]]
                if len(regions_with_pq_1)!=1 or len(regions_with_pq_2)!=1:
                    raise Exception("Each alpha arc has two adjacent regions on two sides.")
                self.region_graph_alpha.add_edge(regions_with_pq_1+regions_with_pq_2+[(a,p),])
        self.region_graph_beta=Graph(len(self.regions),multiedges=True,loops=True)
        for b in self.betas:
            for ind_p,p in enumerate(self.intersections_on_betas[b]):
                q=self.intersections_on_betas[b][ind_p+1-len(self.intersections_on_betas[b])]
                regions_with_pq_1=[R for R in self.regions for ind_foo,foo in enumerate(self.boundary_intersections[R][0::2]) if [foo,self.boundary_intersections[R][2*ind_foo-1]]==[p,q]]
                regions_with_pq_2=[R for R in self.regions for ind_foo,foo in enumerate(self.boundary_intersections[R][0::2]) if [foo,self.boundary_intersections[R][2*ind_foo-1]]==[q,p]]
                if len(regions_with_pq_1)!=1 or len(regions_with_pq_2)!=1:
                    raise Exception("Each beta arc has two adjacent regions on two sides.")
                self.region_graph_beta.add_edge(regions_with_pq_1+regions_with_pq_2+[(b,p),])

        #Now some more error checking. (Definitely not a complete collection.) If diagram is wrong, some error could also be raised by other parts of the program.
        if vector(ZZ,len(self.regions),(len(self.regions))*[1,])*matrix(ZZ,len(self.regions),len(self.intersections),self.point_measures_4)!=vector(ZZ,len(self.intersections),(len(self.intersections))*[4,]):
            raise Exception("Every intersection should have four corners.")
        if len(self.alphas)!=len(self.betas):
            raise Exception("Require same number of alpha and beta circles.")

        #generators is the list of hf-generators, and generator_reps are their representatives. Each generator is represented as a tuple of intersections; the i-th point will be on alpha_i.
        self.generator_reps=[]
        permutations=Arrangements(self.alphas,len(self.alphas)).list()
        for permu in permutations:
            possibilities=[[p for p in self.intersections_on_alphas[a] if self.intersection_incidence[p][1]==permu[a]] for a in self.alphas]
            self.generator_reps+=list(cartesian_product_iterator(possibilities))
        self.generators=range(len(self.generator_reps))#the generator names
        
        self.QSpinC=[] #list whose ith entry is the list of Chern numbers for the spinc structure associated to the ith generator
        for gen in self.generators:
            chern_num=[] #Chern numbers for the spinc structure associated to gen
            gen_vect_abs=vector(ZZ,len(self.intersections))
            for p in self.intersections:
                gen_vect_abs[p]=(self.generator_reps[gen]).count(p)+(self.generator_reps[gen]).count(p)
            for D in self.pd_basis:
                chern_num.append(D.dot_product(self.euler_measures_2_un)/2+(D*self.point_measures_4_un*gen_vect_abs)/4)
            self.QSpinC.append(chern_num)
        #Some more variables that will not be initialized initially, since we need to solve matrix equations.
        self.SpinC="Not yet initialized"#SpinC[i] will be the SpinC structure of the i-th generator (numbered arbitrarily starting at 0).
        self.SpinC_structures="Not yet initialized"#the SpinC structures
        self.genSpinC="Not yet initialized"#genSpinC[i] is a list of generators indices living in the i-th SpinC structure.
        #domains_stored[i][j] stores (m,D) where D is a domain (as a vector over the unpointed regions) from the i-th generator to the j-th generator, and m is its Maslov index. If there are multiple domains (happens when H^1 not 0). If there are none (happens when H_1 not 0, i.e., have multiple SpinC structures), stores None.
        self.domains_stored="Not yet initialized"#domains_stored will be a double list, domains_stored[i][j] will store (m,D) or None; where D is some domain from i to j and m is its Maslov index.
        self.gradings="Not yet initialized"#The absolute gradings of the generators; can be reset manually by set_abs_gr.
        self.pos_differentials=dict()#a dictionary of Maslov index 1 positive domains between pairs of generetors
        self.complex_initialized=False#will be set to true once the chain complex is initialized.
        self.differential_initialized=False#will be set to true once the differential is initialized.
        #The chain complexes will also not be initialized initially.

        #Variables needed to solve for domains.
        k=(self.boundary_mat_un.transpose()).nrows() #the same as the number of intersection points
        l=(self.boundary_mat_un.transpose()).ncols() #the same as the number of unpointed regions
        P, L, U = (self.boundary_mat_un.transpose()).LU() #B=PLU with B=self.boundary_mat_un.transpose()
        A=P*L #kxk invertible matrix
        self.Ainv=A.inverse()
        self.zrlist=[i for i in range(k) if U.row(i) == zero_vector(ZZ,l)] #list of zero rows of U
        self.plist=U.pivots() #pivot columns of U
        self.prlist=U.pivot_rows() #pivot rows of U
        self.nonplist=[i for i in range(l) if i not in self.plist] #non-pivot columns of U
        self.nonprlist=[i for i in range(k) if i not in self.prlist] #non-pivot rows of U
        self.U_red=U.matrix_from_rows_and_columns(self.prlist,self.plist) #extract pivot rows and columns of U
        self.Ainv_test=(self.Ainv).matrix_from_rows(self.zrlist) #extract the rows of Ainv corresponding to the zero rows of U
        self.Ainv_red=(self.Ainv).delete_rows(self.nonprlist) #remove the rows of Ainv corresponding to the nonpivot rows of U
        self.U_red_inv=(self.U_red).inverse() #U_red is an invertible (k-z)x(k-z) matrix where z is the length of nonprlist
        self.solve_mat=self.U_red_inv*self.Ainv_red #(k-z)xk matrix
        self.red_em_2_un=matrix(self.euler_measures_2[:len(self.regions_un)])
        self.red_em_2_un=self.red_em_2_un.delete_columns(self.nonplist)
        self.red_em_2_un_vect=vector(self.red_em_2_un)
        self.red_pm_4_un=matrix(self.point_measures_4[:len(self.regions_un)])
        self.red_pm_4_un=self.red_pm_4_un.delete_rows(self.nonplist)
        
    def cycle_reps(self,g):
        #Determines permutation type of a generator in the Heegaard diagram of an open book decomposition
        gen=self.generator_reps[g]
        cycle_gen=[]
        while len(flatten(cycle_gen))<len(gen):
            start_p=next(p for p in gen if p not in flatten(cycle_gen))
            new_cycle=[]
            curr_p=start_p
            while curr_p!=start_p or new_cycle==[]:
                new_cycle.append(curr_p)
                found_next_point=False
                next_p=next(p for p in gen if self.intersection_incidence[curr_p][1]==self.intersection_incidence[p][0])
                found_next_point=True
                curr_p=next_p
            if not found_next_point:
                curr_p=start_p
            cycle_gen.append(new_cycle)
        return cycle_gen
    
    def find_domain(self,initial,final):
        #Finds a domain D from the initial generator to the final generator; returns (Maslov_index(D),D). Called by sort_generators.
        target_vect=zero_vector(len(self.intersections))
        target_vect_abs=zero_vector(len(self.intersections))
        for p in list(self.generator_reps[final]):
            target_vect[p]+=1
            target_vect_abs[p]+=1
        for p in list(self.generator_reps[initial]):
            target_vect[p]-=1
            target_vect_abs[p]+=1
        if self.Ainv_test*target_vect==zero_vector(len(self.zrlist)):#checks if there is an answer over the rationals
            solution=self.solve_mat*target_vect
            if [i for i in range(len(solution)) if not (solution[i]).is_integer()]==[]:
                maslov_4=2*solution.dot_product(self.red_em_2_un_vect)+solution*self.red_pm_4_un*target_vect_abs
                if maslov_4%4!=0:
                    return None
                else:
                    maslov=maslov_4/4
                    solution_list=list(solution)
                    for i in range(len(self.nonplist)):
                        solution_list.insert(self.nonplist[i],0)
                    domain=vector(ZZ,len(self.regions_un),solution_list)
                    return (maslov,domain)
            else:
                return None
        else:
            return None
    
    def sort_generators(self):
        #Sorts generators into SpinC and Maslov gradings while recording domains between pairs of generators.
        if self.complex_initialized:
            return True
        
        self.domains_stored=[]
        self.SpinC=[] #list whose ith entry is the spinc structure associated to the ith generator as specified by its Chern numbers and its torsion component
        for gen in self.generators:
            if gen==min([foo for foo in self.generators if self.QSpinC[foo]==self.QSpinC[gen]]):
                (self.SpinC).append((self.QSpinC[gen],0))
            else:
                (self.SpinC).append((self.QSpinC[gen],None))
        for initial in self.generators:
            self.domains_stored.append([])
            same_QSpinC=[gen for gen in self.generators if self.QSpinC[gen]==self.QSpinC[initial]] #generators lying in spinc strutures with the same Chern numbers as that of the initial generator
            for final in self.generators:
                if final not in same_QSpinC:
                    (self.domains_stored[initial]).append(None)
                else:
                    S=self.SpinC[final]
                    if final<initial:
                        if self.domains_stored[final][initial]!=None:
                            (m,D)=self.domains_stored[final][initial]
                            (self.domains_stored[initial]).append((-m,-D))
                        else:
                            (self.domains_stored[initial]).append(None)
                    elif final==initial:
                        (self.domains_stored[initial]).append((0,vector(ZZ,len(self.regions_un))))
                    else:
                        if S[1]==None:
                            if self.find_domain(initial,final)!=None: #checks if the initial and final generators are in the same spinc structure
                                (self.domains_stored[initial]).append(self.find_domain(initial,final))
                                self.SpinC[final]=self.SpinC[initial]
                            elif same_QSpinC.index(final)==same_QSpinC.index(initial)+1: #if the final generator is the first generator in some same spinc structure with the same Chern numbers that comes after the initial generator and is not in the same spinc structure
                                self.SpinC[final]=(self.SpinC[final][0],max([self.SpinC[gen][1] for gen in same_QSpinC])+1)
                                (self.domains_stored[initial]).append(None)
                            else:
                                (self.domains_stored[initial]).append(None)
                        elif S==self.SpinC[initial]: #if the initial and final generators are in the same spinc structure
                            first=next(gen for gen in self.generators if self.SpinC[gen]==S)
                            (m1,D1)=self.domains_stored[first][initial]
                            (m2,D2)=self.domains_stored[first][final]
                            (self.domains_stored[initial]).append((m2-m1,D2-D1))
                        else:
                            (self.domains_stored[initial]).append(None)

        self.SpinC_structures=[] #list of spinc structures with non-empty sets of generators
        for S in self.SpinC: 
            if S not in self.SpinC_structures:
                (self.SpinC_structures).append(S)
        self.div=dict() #list of integers whose ith entry is the divisibility of the ith spinc structure
        for S in self.SpinC_structures:
            ind_S=(self.SpinC_structures).index(S)
            if S[0]==[]:
                self.div[ind_S]=0
            else:
                self.div[ind_S]=gcd(S[0])
        self.genSpinC=[[gen for gen in self.generators if self.SpinC[gen]==S] for S in self.SpinC_structures] #list whose ith entry is a list of generators lying in the ith spinc structure
        self.gradings=[]
        for S in self.SpinC_structures:
            ind_S=(self.SpinC_structures).index(S)
            base_gen=self.genSpinC[ind_S][0]
            if S[0]==[0]*self.b1: #if S is a torsion spinc structure
                abs_gr=dict()
                abs_gr[base_gen]=0
                for gen in self.genSpinC[ind_S]:
                    abs_gr[gen]=abs_gr[base_gen]+self.domains_stored[gen][base_gen][0]
                (self.gradings).append(abs_gr)
            else:
                rel_gr=dict()
                rel_gr[base_gen]=0
                for gen in self.genSpinC[ind_S]:
                    rel_gr[gen]=(rel_gr[base_gen]+self.domains_stored[gen][base_gen][0])%self.div[ind_S]
                (self.gradings).append(rel_gr)
                
        self.complex_initialized=True   

    def find_pos_domains(self,initial,final):
        #Finds all positive domains from the initial generator to the final generator; returns a list of (Maslov_index(D),D). Called by find_diffs_from and find_diffs_to.
        pos_domains=[]
        if self.find_domain(initial,final)!=None:
            (m,D)=self.find_domain(initial,final)
            if self.b1==0 and [i for i in self.regions_un if D[i]<0]==[]:
                pos_domains=[self.find_domain(initial,final)]
            else:
                domain=self.find_domain(initial,final)[1]
                reduced_domain=[]
                for region in self.adjacent_region_list:
                    reduced_domain.append(domain[region])
                target_abs=vector(ZZ,len(self.intersections))
                for p in self.intersections:
                    target_abs[p]=(self.generator_reps[final]).count(p)+(self.generator_reps[initial]).count(p)
                multiplicities=[]
                multiplicities.append([0]*len(self.adjacent_region_list))
                indicator_alphas=[]
                for i in self.alphas:
                    if self.contact_intersections[i] in self.generator_reps[final] and self.contact_intersections[i] not in self.generator_reps[initial]:
                        indicator_alphas.append(i)
                if indicator_alphas!=[]:
                    uniq_adjacent_region_list=list(self.adjacent_region_list)
                    flat_adjacent_region_list=flatten(self.region_pairs)
                    truncated_indicator_alphas=[]
                    for a in indicator_alphas:
                        region_a_0=self.region_pairs[a][0]
                        region_a_1=self.region_pairs[a][1]
                        if flat_adjacent_region_list.count(region_a_0)==1 and flat_adjacent_region_list.count(region_a_1)==1:
                            truncated_indicator_alphas.append(a)
                            uniq_adjacent_region_list.remove(region_a_0)
                            uniq_adjacent_region_list.remove(region_a_1)
                        elif flat_adjacent_region_list.count(region_a_0)!=1 and uniq_adjacent_region_list.count(region_a_0)==1:
                            truncated_indicator_alphas.append(a)
                            uniq_adjacent_region_list.remove(region_a_0)
                        elif flat_adjacent_region_list.count(region_a_1)!=1 and uniq_adjacent_region_list.count(region_a_1)==1:
                            truncated_indicator_alphas.append(a)
                            uniq_adjacent_region_list.remove(region_a_1)
                    mult_at_contact=[[1,0],[0,1]]
                    choices=sorted(cartesian_product([[0,1] for i in range(len(truncated_indicator_alphas))]))
                    for choice in choices:
                        mult=[0]*len(self.adjacent_region_list)
                        for a in indicator_alphas:
                            region_a_0=self.region_pairs[a][0]
                            region_a_1=self.region_pairs[a][1]
                            if flat_adjacent_region_list.count(region_a_0)==1 and flat_adjacent_region_list.count(region_a_1)==1:
                                a_ind=truncated_indicator_alphas.index(a)
                                mult[(self.adjacent_region_list).index(region_a_0)]=mult_at_contact[choice[a_ind]][0]
                                mult[(self.adjacent_region_list).index(region_a_1)]=mult_at_contact[choice[a_ind]][1]
                            elif flat_adjacent_region_list.count(region_a_0)!=1:
                                if a in truncated_indicator_alphas:
                                    a_ind=truncated_indicator_alphas.index(a)
                                    mult[(self.adjacent_region_list).index(region_a_0)]=mult_at_contact[choice[a_ind]][0]
                                mult[(self.adjacent_region_list).index(region_a_1)]=1-mult[(self.adjacent_region_list).index(region_a_0)]
                            elif flat_adjacent_region_list.count(region_a_1)!=1:
                                if a in truncated_indicator_alphas:
                                    a_ind=truncated_indicator_alphas.index(a)
                                    mult[(self.adjacent_region_list).index(region_a_1)]=mult_at_contact[choice[a_ind]][1]
                                mult[(self.adjacent_region_list).index(region_a_0)]=1-mult[(self.adjacent_region_list).index(region_a_1)]
                            else:
                                pass
                        if mult not in multiplicities:
                            multiplicities.append(mult)       
                m=len(self.adjacent_region_list)
                K=len(multiplicities)
                mult_matrix=Matrix(ZZ,m,K)
                for i in range(m):
                    for j in range(K):
                        mult_matrix[i,j]=multiplicities[j][i]
                if self.pd_basis!=[]:
                    domain_matrix=Matrix(ZZ,len(self.adjacent_region_list),len(multiplicities))
                    for i in range(len(self.adjacent_region_list)):
                        for j in range(len(multiplicities)):
                            domain_matrix[i,j]=reduced_domain[i]
                    delta_matrix=mult_matrix-domain_matrix
                    for k in range(len(multiplicities)):
                        delta=delta_matrix.column(k)
                        try:    
                            coef=(self.indicator_matrix).solve_right(delta)
                        except:
                            pass
                        else:
                            pD=[0]*len(self.regions_un)
                            for j in range(len(self.regions_un)):
                                for i in range(len(self.pd_basis)):
                                    pD[j]+=coef[i]*self.pd_basis[i][j]
                            new_domain=[domain[j]+pD[j] for j in range(len(self.regions_un))]
                            new_domain=vector(ZZ,new_domain)
                            if [a for a in new_domain if a<0]==[]:
                                maslov_4=2*new_domain.dot_product(self.euler_measures_2_un)+new_domain.dot_product(self.point_measures_4_un*target_abs)
                                maslov=maslov_4/4
                                pos_domains.append((maslov,new_domain))
        else:
            pos_domains=[]
        return pos_domains

    def domain_type(self,initial,final,domain):
        #Computes J_+, the Euler characteristic, the number of boundary components, and the genus of a Lipshitz curve associated to a given domain.
        target_abs=vector(ZZ,len(self.intersections))
        for p in self.intersections:
            target_abs[p]=(self.generator_reps[final]).count(p)+(self.generator_reps[initial]).count(p)
        num_trivial=len([a for a in self.generator_reps[initial] if a in self.generator_reps[final]]) #Number of trivial disks
        alphas=list(self.alphas)
        bdy_comps=[]
        while alphas!=[]:
            bdy_comp=[]
            first_alpha=min(alphas)
            current_alpha=first_alpha
            reached_end_of_loop=False
            while not reached_end_of_loop:    
                next_alpha=next(a for a in self.alphas if self.intersection_incidence[self.generator_reps[final][current_alpha]][1]==self.intersection_incidence[self.generator_reps[initial][a]][1])
                bdy_comp.append(next_alpha)           
                alphas.remove(next_alpha)
                current_alpha=next_alpha
                if current_alpha==first_alpha:
                    reached_end_of_loop=True
                    bdy_comps.append(bdy_comp)
        maslov_4=2*domain.dot_product(self.euler_measures_2_un)+domain.dot_product(self.point_measures_4_un*target_abs)
        maslov=maslov_4/4
        euler_measure_2=domain.dot_product(self.euler_measures_2_un)
        iota=maslov-euler_measure_2
        jplus=iota+len(self.cycle_reps(initial))-len(self.cycle_reps(final))
        euler_char=len(self.alphas)-num_trivial-iota #Euler characteristic of the Lipshitz curve
        boundary_comps=len(bdy_comps)-num_trivial
        genus=(2-euler_char-boundary_comps)/2 #Genus of the Lipshitz curve
        
        return(jplus,euler_char,boundary_comps,genus)
    
    def can_contribute(self,initial,final):
        #Checks if there is a Maslov index 1 positive domain from the initial generator to the final generator.
        if not self.complex_initialized:
            self.sort_generators()
        if initial!=final and self.find_pos_domains(initial,final)!=[]:
            self.pos_differentials[(initial,final)]=[]
            for pos_domain in self.find_pos_domains(initial,final):
                if pos_domain[0]==1:
                    self.pos_differentials[(initial,final)].append(pos_domain[1])
            if self.pos_differentials[(initial,final)]!=[]:
                return True
            else:
                return False
        else:
            return False
    
    def find_pos_diffs_from(self,initial):
        #Finds Maslov index 1 positive domains D starting at a given initial generator; returns a list of final generators.
        if not self.complex_initialized:
            self.sort_generators()
        S=self.SpinC[initial]
        same_SpinC=self.genSpinC[(self.SpinC_structures).index(S)] #generators lying in the same spinc struture
        diffs_from=[]

        for final in same_SpinC:
            if final==initial:
                pass
            elif self.pd_basis==[]:
                try:
                    (m,D)=tuple(self.find_domain(initial,final))
                    if m==1 and [a for a in D if a<0]==[]:
                        diffs_from.append(final)
                except:
                    pass
            else:
                for k in range(len(self.find_pos_domains(initial,final))):
                    try:
                        (m,D)=tuple(self.find_pos_domains(initial,final)[k])
                        if m==1 and [a for a in D if a<0]==[]:
                            diffs_from.append((D,final))
                    except:
                        pass
        return diffs_from
    
    def find_pos_diffs_to(self,final):
        #Finds Maslov index 1 positive domains D ending in a given final generator; returns a list of initial generators.
        if not self.complex_initialized:
            self.sort_generators()
        S=self.SpinC[final]
        same_SpinC=self.genSpinC[(self.SpinC_structures).index(S)] #generators lying in the same spinc struture
        diffs_to=[]

        for initial in same_SpinC:
            if initial==final:
                pass
            elif self.pd_basis==[]:
                try:
                    (m,D)=tuple(self.find_domain(initial,final))
                    if m==1 and [a for a in D if a<0]==[]:
                        diffs_to.append(initial)
                except:
                    pass
            else:
                for k in range(len(self.find_pos_domains(initial,final))):
                    try:
                        (m,D)=tuple(self.find_pos_domains(initial,final)[k])
                        if m==1 and [a for a in D if a<0]==[]:
                            diffs_to.append((D,initial))
                    except:
                        pass
        return diffs_to

    def print_differentials(self):
        #Returns and records all Maslov index 1 positive domains along with their topological types.
        if not self.complex_initialized:
            self.sort_generators()
        ft=open(str(self.name)+'_possible_differentials.txt','a')
        for S in self.SpinC_structures:
            ind_S=(self.SpinC_structures).index(S)
            if len(self.genSpinC[ind_S])>1:
                ft.write("In SpinC structure "+str(repr(S))+'\n')
                for initial in self.genSpinC[ind_S]:
                    for final in [gen for gen in self.genSpinC[ind_S] if self.can_contribute(initial,gen) and gen!=initial]:
                        ft.write(str(repr(initial))+" in (relative) grading "+str(repr(self.gradings[ind_S][initial]))+" has possible differentials to "+str(repr(final))+" with domain types (Euler char, genus, num of boundaries): "+str(repr([(self.domain_type(initial,final,domain)[1],self.domain_type(initial,final,domain)[3], self.domain_type(initial,final,domain)[2]) for domain in self.pos_differentials[(initial,final)]]))+'\n')
        ft.close()
        self.differential_initialized=True 
        return None
            
    def plot_complex(self,S,output_file=None):
        #Plots the chain complex in the spinc structure S.
        if not self.differential_initialized:
            self.print_differentials()
        diffs=[]
        colors={'1': "blue", '(?)': "red"}
        heights=self.gradings
        for initial in self.genSpinC[S]:
            for final in [gen for gen in self.genSpinC[S] if self.can_contribute(initial,gen) and gen!=initial]:
                if self.pos_differentials[(initial,final)]!=[]:
                    for domain in self.pos_differentials[(initial,final)]:
                        if self.domain_type(initial,final,domain)[1]==1 and self.domain_type(initial,final,domain)[2]==1:
                            diffs.append((initial,final,'1'))
                        else:
                            diffs.append((initial,final,'(?)'))
        G = DiGraph(diffs, multiedges=True) 
        Gplot = plot(G, edge_labels=false, edge_colors=G._color_by_label(colors), layout='acyclic', edge_labels_background= 'transparent', figsize=[40, 20], vertex_size=100)
        if output_file==None:
            Gplot.show() 
        else:
            Gplot.save(output_file)
    
    #Use the following functions to analyze the canonical (Honda-Kazez-Matic) generator
    def find_zero_diffs_from(self,initial):
        #Finds J_+=0, Maslov index 1 positive domains starting at a given initial generator; returns a list of final generators.
        if not self.complex_initialized:
            self.sort_generators()
        S=self.SpinC[initial]
        same_SpinC=self.genSpinC[(self.SpinC_structures).index(S)] #generators lying in the same spinc struture
        zero_diffs_from=[]

        for final in same_SpinC:
            if final==initial:
                pass
            elif self.pd_basis==[]:
                try:
                    (m,D)=tuple(self.find_domain(initial,final))
                    if m==1 and [a for a in D if a<0]==[] and self.domain_type(initial,final,D)[0]==0:
                        zero_diffs_from.append(final)
                except:
                    pass
            else:
                for k in range(len(self.find_pos_domains(initial,final))):
                    try:
                        (m,D)=tuple(self.find_pos_domains(initial,final)[k])
                        if m==1 and [a for a in D if a<0]==[] and self.domain_type(initial,final,D)[0]==0:
                            zero_diffs_from.append(final)
                    except:
                        pass
        return zero_diffs_from
    
    def find_zero_diffs_to(self,final):
        #Finds J_+=0, Maslov index 1 positive domains ending in a given final generator; returns a list of initial generators.
        if not self.complex_initialized:
            self.sort_generators()
        S=self.SpinC[final]
        same_SpinC=self.genSpinC[(self.SpinC_structures).index(S)] #generators lying in the same spinc struture
        zero_diffs_to=[]

        for initial in same_SpinC:
            if initial==final:
                pass
            elif self.pd_basis==[]:
                try:
                    (m,D)=tuple(self.find_domain(initial,final))
                    if m==1 and [a for a in D if a<0]==[] and self.domain_type(initial,final,D)[0]==0:
                        zero_diffs_to.append(initial)
                except:
                    pass
            else:
                for k in range(len(self.find_pos_domains(initial,final))):
                    try:
                        (m,D)=tuple(self.find_pos_domains(initial,final)[k])
                        if m==1 and [a for a in D if a<0]==[] and self.domain_type(initial,final,D)[0]==0:
                            zero_diffs_to.append(initial)
                    except:
                        pass
        return zero_diffs_to
    
    def zero_crawl(self):
        #Builds the smallest subcomplex in the canoncial spinc structure containing the contact generator and only those differentials with J_+=0.
        self.zero_diffs=[]
        self.cx00=[0]
        self.cx10=[]
        new_output=[0]
        control=True
        while control==True:
            new_input=[]
            for gen2 in new_output:
                for gen1 in self.find_zero_diffs_to(gen2):
                    if (gen1,gen2) not in self.zero_diffs:
                        (self.zero_diffs).append((gen1,gen2))
                    if gen1 not in self.cx10:
                        (self.cx10).append(gen1)
                        new_input.append(gen1)
            if new_input!=[]:
                new_output=[]
                for gen1 in new_input:
                    for gen2 in self.find_zero_diffs_from(gen1):
                        if (gen1,gen2) not in self.zero_diffs:
                            (self.zero_diffs).append((gen1,gen2))
                        if gen2 not in self.cx00:
                            (self.cx00).append(gen2)
                            new_output.append(gen2)
                if new_output==[]:
                    control=False
            else: 
                control=False
    
    def check_zero_order(self):
        #Determines if the EH class is zero via J_+=0 differentials
        self.cx=self.genSpinC[0]
        self.cx0=[]
        self.cx1=[]
        for gen in self.cx:
            if self.gradings[0][gen]==0:
                (self.cx0).append(gen)
            elif self.gradings[0][gen]==1:
                (self.cx1).append(gen)
        B=Matrix(GF(2),len(self.cx0),len(self.cx1))
        for final in self.cx0:
            for initial in self.find_zero_diffs_to(final):
                total=self.find_zero_diffs_to(final).count(initial)
                i=(self.cx0).index(final)
                j=(self.cx1).index(initial)
                if total%2==0:
                    B[i,j]=0
                else:
                    B[i,j]=1
        m=B.nrows()
        n=B.ncols()
        Y=matrix(GF(2),m,1)
        for i in range(m):
            if i==0:
                Y[i,0]=1
            else:
                Y[i,0]=0
        try:
            X=B.solve_right(Y)
            self.killer_chain=[self.cx1[i] for i in range(n) if X[i]!=0]
            print("Order is zero via "+str(self.killer_chain))
        except:
            print("Order is non-zero")