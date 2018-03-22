from arcreactor.protobufs.graph_pb2 import Graph
from arcreactor.protobufs.kinetics_pb2 import SystemKinetics as Kinetics
import arcreactor, copy, asyncio
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage

class Reactors:
    
    PFR = 'pfr'
    CSTR = 'cstr'
    
    def __init__(self):
        self.graph = Graph()
        source = self.graph.nodes[0]
        source.id = 0
        source.label = 'source'
        self.node_ids = 1
        self.edge_ids = 0
        self._update_nxgraph()
        self._reset()
    
    def _reset(self):
        self.system = arcreactor.simulation.Simulation(0)
        self.state = Kinetics()
        
    def _step(self, dt=1):
        tmp_graph = copy.copy(self.graph)
        tmp_state = copy.copy(self.state)
        self.system.graph_time += int(dt)
        loop = asyncio.get_event_loop()
        self.state = loop.run_until_complete(self.system.calculate(tmp_state, tmp_graph))
            
    def add_reactor(self, reactor_type, temperature = 300):
        assert reactor_type == 'cstr' or reactor_type == 'pfr'
        node = self.graph.nodes[self.node_ids]
        node.id = self.node_ids
        node.label = reactor_type
        node.weight.append(temperature)
        self.node_ids += 1
        node.delete = False
        self._update_nxgraph()
        return node.id
    
    def connect(self, inlet, outlet):
        assert len(self.graph.nodes) > inlet and len(self.graph.nodes) > outlet
        assert not self.graph.nodes[inlet].delete and not self.graph.nodes[outlet].delete
        edge = self.graph.edges[self.edge_ids]
        edge.idA = inlet
        edge.idB = outlet
        edge.labelA = self.graph.nodes[edge.idA].label
        edge.labelB = self.graph.nodes[edge.idB].label
        self.edge_ids += 1
        self._update_nxgraph()
    
    @property
    def nxgraph(self):
        if self._nxgraph is None:
            self._update_nxgraph()
        return self._nxgraph
    
    def _update_nxgraph(self):
        G = nx.DiGraph()
        for i in range(len(self.graph.nodes)):
            n = self.graph.nodes[i]
            if not n.delete:
                G.add_node(n.id, label=n.label)
        for i in range(len(self.graph.edges)):
            e = self.graph.edges[i]
            if not self.graph.nodes[e.idA].delete and not self.graph.nodes[e.idB].delete:
                G.add_edge(e.idA, e.idB)
        self._nxgraph = G


    def _layout_graph(self):
        layout = nx.nx_agraph.graphviz_layout(self.nxgraph, prog='dot')        
        return layout
    
    def _plot_graph(self, layout, ax, radius=100, offset=30):
        #radius is in an unknown unit system. 
        #networkx defaults to 300
        colors = []
        labels = {}
        for i in range(len(self.graph.nodes)):
            n = self.graph.nodes[i]
            if n.delete:
                continue
            if n.label == 'source':
                colors.append('gray')
                labels[n.id] = 'Source'
            elif n.label == 'pfr':
                colors.append('white')
                labels[n.id] = 'PFR'
            else:
                colors.append('black')
                labels[n.id] = 'CSTR'            
        nx.draw_networkx_nodes(self.nxgraph, pos=layout, ax=ax, node_size=radius, 
                                     node_color=colors)
        nx.draw_networkx_edges(self.nxgraph, pos=layout, ax=ax)
        #pos mod is to make labels appear below center
        nx.draw_networkx_labels(self.nxgraph, pos={ k:(x[0], x[1] - offset) for k,x in layout.items()}, labels=labels)
    def _plot_fracs(self, layout, ax, radius=50):
        ax.axis('equal')
        #set-up colors
        N = len(self.state.kinetics[0].mole_fraction)
        cmap = plt.get_cmap('Accent')
        colors = [cmap(i / (N - 1)) for i in range(N)]
        #plot source
        ax.add_artist(mpl.patches.Circle(layout[0], radius, facecolor='darkgray'))
        for i in range(len(self.state.kinetics)):
            #check for unconnected reactors
            if sum(self.state.kinetics[i].mole_fraction) < 0.1:
                ax.add_artist(mpl.patches.Circle(layout[self.state.kinetics[i].id], radius, edgecolor='darkgray', facecolor='white'))
            else:
                #why do things not add up to 1 sometimes?
                plt.pie([x / sum(self.state.kinetics[i].mole_fraction) for x in self.state.kinetics[i].mole_fraction],
                    radius=radius, 
                    center=layout[self.state.kinetics[i].id], 
                    frame=True,
                    colors = colors,
                    autopct = lambda x : '{:.0f}%'.format(x),
                    pctdistance  = 0.8)
    def plot_reactors(self, time = 0, fig = None):
        if fig is None:
            fig, ax = plt.subplots()
        else:
            ax = fig.gca()
        self._reset()
        for i in range(time):
            self._step()
        layout = self._layout_graph()
        self._plot_fracs(layout, ax, radius=25)
        self._plot_graph(layout, ax)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.legend([r'H$_2$', r'CO$_2$', r'O$_2$', r'H$_2$O'])
        ax.set_title('$t = {:.2f}$'.format(self.system.graph_time))
        plt.show()
    
    def animate_reactors(self, time, fps = 30):
        fig = plt.figure()
        ax = fig.gca()
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        
        self._reset()
        
        layout = self._layout_graph()
        self._step(0)
        
        duration = time / fps + 3
        
        def draw(t):
            if t > 1 and duration - t > 2:
                #last 2 seconds and first sceond are filler
                self._step() 
            ax.clear()
            self._plot_fracs(layout, ax, radius=25)
            self._plot_graph(layout, ax)
            ax.legend([r'H$_2$', r'CO$_2$', r'O$_2$', r'H$_2$O'])
            ax.set_title('$t = {:.2f}$'.format(self.system.graph_time))
            return mplfig_to_npimage(fig)
        animation = VideoClip(draw, duration=duration)
        return animation