use std::{io, iter, fmt, str};
use std::usize;
use std::num::{ParseIntError, ParseFloatError};
use std::io::BufRead;

// svg things
const SVG_LEAVE_TOP: f32 = 27.;
const SVG_BLINDER_PADDING_TOP: f32 = 27.;
const SVG_BLINDER_PADDING_BOTTOM: f32 = 20.;
const GRAPH_WIDTH: f32 = 450.;
const LEVEL_HEIGHT: f32 = 50.;
const GROUP_COLORS: &[&str] = &[
    "dodgerblue",
    "darkorange",
    "mediumvioletred",
];

/// ROUTINE LABEL
fn routine_label(d: &Graph) -> Vec<usize> {
    // %#%- 1.
    let (
        // what nodes are in [level]?
        nodes_in_level,
        // in what level is [node]?
        level,
    ) = determine_node_level(d);

    // a linear labelling. label=lin_label[node]
    let mut lin_label = vec![usize::MAX; d.nodes.len()];
    let mut rev_lin_label = Vec::with_capacity(d.nodes.len());
    for level in &nodes_in_level {
        for &node in level {
            let label = rev_lin_label.len();
            rev_lin_label.push(node);
            lin_label[node] = label;
        }
    }

    // how many levels are there?
    let number_of_levels = nodes_in_level.len();

    svg_start("LEVEL OF EACH NODE");
    println!("<!-- levels -->");
    for (node, label) in level.iter().enumerate() {
        let node = &d.nodes[node];
        println!(r#"<text font-style="italic" text-anchor="start" x="{}" y="{}">{}</text>"#,
            node.position[0] + 8., node.position[1], label+1);
    }
    svg_end(&d);

    // %#%- 2. (initialization)
    // Referred to in the paper as EDGES.
    let mut edges_leaving_level = vec![vec![]; number_of_levels];
    // what node has [label]?
    let mut node_by_label = vec![usize::MAX; d.nodes.len()];
    // what label does [node] have?
    let mut label = vec![usize::MAX; d.nodes.len()];
    // labels start at 0
    let mut j = 0;

    // %#%- 3 and 4.
    for i in 0..number_of_levels {
        let j_old = j;
        routine_partition(i, &mut edges_leaving_level, &mut j, &mut node_by_label, &mut label, &nodes_in_level, d, &lin_label, &rev_lin_label);
        // %#%- 3.1.
        // for each node y at level i, in order of increasing label
        for y in (j_old..j).map(|l| node_by_label[l]) {
            for &x in &d.reverse_edges[y] {
                edges_leaving_level[level[x]].push((x, y));
            }
        }
    }

    label
}

fn routine_partition(i: usize, edges_leaving_level: &mut Vec<Vec<(usize, usize)>>, j: &mut usize, node_by_label: &mut Vec<usize>, label: &mut Vec<usize>, nodes_in_level: &Vec<Vec<usize>>, graph: &Graph, _lin_label: &Vec<usize>, rev_lin_label: &Vec<usize>){
    // NOTE: All edges in EDGES have a source node in level i
    // NOTE: For each pair of edges, let y and z be the destination nodes.
    //     y will be higher to the top of the stack than z if L(y) > L(z).
    // - - - - - -
    struct Succ{
        nodes: Vec<usize>,
        low_label: usize,
        temp_list: Vec<usize>,
    }

    let num_nodes_in_level_i = nodes_in_level[i].len();
    // used to remove elements from Succ::nodes in constant time.
    let mut position_in_group_deck = vec![usize::MAX; num_nodes_in_level_i];

    // %#%- 1.
    let mut all_succ = vec![Succ{low_label: *j, temp_list: Vec::new(), nodes: Vec::new() }];
    let mut succ = vec![usize::MAX; num_nodes_in_level_i];
    for &node in &nodes_in_level[i] {
        let w = 0;
        position_in_group_deck[rev_lin_label[node]-*j] = all_succ[w].nodes.len();
        all_succ[w].nodes.push(node);
        succ[rev_lin_label[node]-*j] = w;
    }

    // %#%- 2.
    let elli = &mut edges_leaving_level[i];
    while let Some(&(_, y)) = elli.last() {
        // used for 2.1.2.
        let mut touched_temp_lists = Vec::new();

        svg_start(&format!(
            "{}{} LEVEL WIP GROUPINGS",
            i+1,
            if i==0 { "ST"} else if i == 1 { "ND" } else if i == 2 { "RD" } else {"TH"}
        ) );

        println!("<!-- LINES TO Y -->");
        let to = &graph.nodes[y];
        println!(r#"<circle cx="{}" cy="{}" r="24" fill="none" stroke=springgreen stroke-width="4" opacity=0.8 />"#,
            to.position[0], to.position[1] - 5.
        );

        println!("<!-- WIP GROUPINGS -->");
        let mut discount = 0;
        for (s,succ) in all_succ.iter().enumerate() {
            let sd = s - discount;
            if succ.nodes.len() == 0 {
                discount += 1;
                continue;
            }
            for &node in &succ.nodes {
                let node = &graph.nodes[node];
                println!(r#"<circle cx="{}" cy="{}" r="5" fill="{}"/>"#,
                    node.position[0], node.position[1] - 5., GROUP_COLORS[sd]);
                println!(r#"<text text-anchor="start" x="{}" y="{}" font-size="small">{}=LL</text>"#,
                    node.position[0] + 10., node.position[1] - 3., succ.low_label + 1
                );
            }
        }

        // %#%- 2.1.
        // %#%- 2.1.1.
        while let Some(&(x, y2)) = elli.last() {
            // %#%- 2.1.1.1.
            if y != y2 { break; }
            let w = succ[rev_lin_label[x]-*j];
            let posdeck = &mut position_in_group_deck[rev_lin_label[x]-*j];

            let &displaced = all_succ[w].nodes.last().unwrap();
            all_succ[w].nodes[*posdeck] = displaced;
            all_succ[w].nodes.pop();
            *posdeck = all_succ[w].temp_list.len();
            position_in_group_deck[rev_lin_label[displaced]-*j] = *posdeck;
            all_succ[w].temp_list.push(x);
            touched_temp_lists.push(w);
            elli.pop();
            
            let from = &graph.nodes[x];
            println!(r#"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="springgreen" stroke-width="5" opacity=0.8 />"#,
                from.position[0], from.position[1] - 5.,
                to.position[0],   to.position[1]   - 5.);
        }

        println!("<!-- BLINDERS -->");
        if graph.nodes[nodes_in_level[i][0]].position[1] - SVG_LEAVE_TOP - SVG_BLINDER_PADDING_TOP > 0. {
            println!( r#"<rect x="0" y="{}" width="{}" height="{}" fill="url(#hazard)" opacity="0.3" />"#,
                SVG_LEAVE_TOP, GRAPH_WIDTH,
                graph.nodes[nodes_in_level[i][0]].position[1] - SVG_LEAVE_TOP - SVG_BLINDER_PADDING_TOP
            );
        }
        println!( r#"<rect x="0" y="{}" width="{}" height="100%" fill="black" opacity="0.15" />"#,
            graph.nodes[nodes_in_level[i][0]].position[1] + SVG_BLINDER_PADDING_BOTTOM,
            GRAPH_WIDTH
        );
        println!("<!-- ALREADY ASSIGNED LABELS -->");
        for lidx in 0..i {
            for &node in &nodes_in_level[lidx] {
                let l = label[node];
                let node = &graph.nodes[node];
                println!(r#"<text font-style="italic" text-anchor="start" x="{}" y="{}">{}</text>"#,
                    node.position[0] + 8., node.position[1], l+1);
            }
        }
        svg_end(&graph);

        // %#%- 2.1.2.
        for w in touched_temp_lists {
            let v = all_succ.len();
            let nodes = std::mem::take(&mut all_succ[w].temp_list);
            for x in &nodes {
                succ[rev_lin_label[*x]-*j] = v;
            }
            all_succ.push(Succ{
                nodes,
                low_label: all_succ[w].low_label + all_succ[w].nodes.len(),
                temp_list: Vec::new(),
            });
        }
    }

    svg_start(&format!(
        "{}{} LEVEL FINAL GROUPINGS",
        i+1,
        if i==0 { "ST"} else if i == 1 { "ND" } else if i == 2 { "RD" } else {"TH"}
    ) );
    println!("<!-- FINAL GROUPINGS -->");
    let mut discount = 0;
    for (s,succ) in all_succ.iter().enumerate() {
        let sd = s - discount;
        if succ.nodes.len() == 0 {
            discount += 1;
            continue;
        }
        for &node in &succ.nodes {
            let node = &graph.nodes[node];
            println!(r#"<circle cx="{}" cy="{}" r="5" fill="{}"/>"#,
                node.position[0], node.position[1] - 5., GROUP_COLORS[sd]);
            println!(r#"<text font-weight="bold" text-anchor="start" x="{}" y="{}" font-size="small">{}=LL</text>"#,
                node.position[0] + 10., node.position[1] - 3., succ.low_label + 1
            );
        }
    }
    println!("<!-- BLINDERS -->");
    if graph.nodes[nodes_in_level[i][0]].position[1] - SVG_LEAVE_TOP - SVG_BLINDER_PADDING_TOP > 0. {
        println!( r#"<rect x="0" y="{}" width="{}" height="{}" fill="url(#hazard)" opacity="0.3" />"#,
            SVG_LEAVE_TOP, GRAPH_WIDTH,
            graph.nodes[nodes_in_level[i][0]].position[1] - SVG_LEAVE_TOP - SVG_BLINDER_PADDING_TOP
        );
    }
    println!( r#"<rect x="0" y="{}" width="{}" height="100%" fill="black" opacity="0.3" />"#,
        graph.nodes[nodes_in_level[i][0]].position[1] + SVG_BLINDER_PADDING_BOTTOM,
        GRAPH_WIDTH
    );
    println!("<!-- ALREADY ASSIGNED LABELS -->");
    for lidx in 0..i {
        for &node in &nodes_in_level[lidx] {
            let l = label[node];
            let node = &graph.nodes[node];
            println!(r#"<text font-style="italic" text-anchor="start" x="{}" y="{}">{}</text>"#,
                node.position[0] + 8., node.position[1], l+1);
        }
    }
    svg_end(&graph);

    // %#%- 3.
    for &x in &nodes_in_level[i] {
        let s = &mut all_succ[succ[rev_lin_label[x]-*j]];
        label[x] = s.low_label;
        node_by_label[s.low_label] = x;
        s.low_label += 1;
    }

    *j += num_nodes_in_level_i;
}
/// Information on how to draw this Node
///
/// Parsed from the input file
struct Node {
    name: String,
    position: [f32; 2],
    draw_circle: bool,
}

/// Input to the problem
struct Graph {
    nodes: Vec<Node>,
    /// let to_nodes: Vec<_> = edges[from_node];
    edges: Vec<Vec<usize>>,
    /// let from_nodes: Vec<_> = reverse_edges[to_node];
    reverse_edges: Vec<Vec<usize>>,
}

impl Graph {
    pub fn try_parse(read: &mut dyn BufRead) -> Result<Graph, InputError>
    {
        use InputError::*;

        let line = Self::get_line(read)?.ok_or(ShortRead)?;
        if line != "MLDAG" {
            return Err(ParseError);
        }

        let line = Self::get_line(read)?.ok_or(ShortRead)?;
        let num_nodes = line.parse::<usize>().map_err(|e| CantParseNumber(line.to_owned(), e))?;

        let mut nodes = Vec::with_capacity(num_nodes);
        let mut node_by_name = std::collections::HashMap::new();
        for i in 0..num_nodes {
            let name = Self::get_identifier(read)?.ok_or(ShortRead)?;

            let posx = Self::get_identifier(read)?.ok_or(ShortRead)?;
            let posx = posx.parse::<f32>().map_err(|e| CantParseFloat(posx.to_owned(), e))?;

            let posy = Self::get_identifier(read)?.ok_or(ShortRead)?;
            let posy = posy.parse::<f32>().map_err(|e| CantParseFloat(posy.to_owned(), e))?;

            if node_by_name.insert(name.clone(), i).is_some() {
                return Err(NodesHaveSameName(name));
            }

            nodes.push(Node{
                name,
                position: [posx, posy],
                draw_circle: true,
            });
        }

        let line = Self::get_line(read)?.ok_or(ShortRead)?;
        let num_edges = line.parse::<usize>().map_err(|e| CantParseNumber(line.to_owned(), e))?;

        let mut edges = vec![vec![]; num_nodes];
        let mut reverse_edges = vec![vec![]; num_nodes];
        for _ in 0..num_edges {
            let from = Self::get_identifier(read)?.ok_or(ShortRead)?;
            let from = if let Some(&from) = node_by_name.get(&from) { from }
            else { return Err(EdgeReferencesNonexistantNode(from.to_owned())) };

            let to = Self::get_identifier(read)?.ok_or(ShortRead)?;
            let to = if let Some(&to) = node_by_name.get(&to) { to }
            else { return Err(EdgeReferencesNonexistantNode(to.to_owned())) };

            edges[from].push(to);
            reverse_edges[to].push(from);
        }

        Ok(Graph{
            nodes,
            edges,
            reverse_edges
        })
    }
    /// Private helper function for reading a line from a bufread,
    /// while trimming comments, whitespace, and empty lines
    fn get_line(read: &mut dyn BufRead) -> io::Result<Option<String>>
    {
        // read_line returns 1 for empty lines, as the newlines are included
        // and so, as it says in its docs, 0 is EOF.
        let mut read_buf = String::new();
        while read.read_line(&mut read_buf)? != 0 {
            let mut line = read_buf.as_str();

            // remove comments
            // (everything following two slashes)
            if let Some(non_comment_len) = line.find("//") {
                line = &line[0..non_comment_len];
            }

            // remove whitespace
            line = line.trim();

            // return first non-empty line
            if !line.is_empty() {
                return Ok(Some(line.to_string()));
            }

            // Gotta clear read_buf because read_line appends
            read_buf.clear();
        }
        Ok(None)
    }

    /// Private helper function for reading a single identifier from a bufread.
    ///
    /// Identifiers are space delimited sequences of one or more characters
    fn get_identifier(read: &mut dyn BufRead) -> io::Result<Option<String>>
    {
        let mut identifier = String::new();
        // number of bytes we've consumed
        let mut consumed_amount = 0usize;
        // are we currently throwing out data until we hit the newline character?
        let mut in_comment = false;

        // This loop parses one identifier
        // NOTE: This 'ident: syntax is for a named break, and is not a goto.
        'ident: loop {
            // Consume character(s) from previous go through the loop
            if consumed_amount != 0 {
                read.consume(consumed_amount);
                consumed_amount = 0;
            }

            if in_comment {
                in_comment = false;
                let mut comment = String::new();
                read.read_line(&mut comment)?;

                //  if we have an identifier, it's completed by the whitespace
                if !identifier.is_empty() { break 'ident; }
            }

            // buffer of u8's (not unicode 'char's!)
            let buffer = read.fill_buf()?;
            if buffer.is_empty() { break 'ident; }

            let valid_amount = match str::from_utf8(&buffer) {
                // the whole buffer is valid UTF-8
                Ok(_) => buffer.len(),
                // there's a unicode error at the start of the buffer, so it's a real error.
                Err(e) if e.valid_up_to() == 0 =>
                    // error message akin to BufRead::read_line
                    return Err(io::Error::new(io::ErrorKind::InvalidData, "stream did not contain valid UTF-8")),
                // there's an error, but we read something
                Err(e) => e.valid_up_to()
            };


            //
            // we read something; try to add it to the identifier
            //

            // Unwrap will never panic because the previous from_utf call says these are valid bytes
            let valid = str::from_utf8(buffer.split_at(valid_amount).0).unwrap();

            // Iterate through pairs of characters in the input,
            for (c, n) in valid.chars().zip(
                //  including the last character paired with no character following
                valid.chars().map(|c| Some(c)).chain(iter::once(None))
            ) {
                if c=='/' && n == Some('/') {
                    // Nuke until newline or EOF
                    in_comment = true;
                    consumed_amount += 2; // sizeof("//")
                    continue 'ident;
                } else if c=='/' && n == None && valid_amount != 1 {
                    // Can't consume c without next character

                    // NOTE: if valid_amount is 1, then we're at EOF, or the next character is invalid UTF:
                    //   if it's invalid UTF, we're gonna return an error
                    //   if it's EOF, we can consume c.

                    break;
                } else if c.is_whitespace() { // NOTE: newlines count for is_whitespace
                    // Consume c without appending to identifier,
                    consumed_amount += c.len_utf8();
                    //  if we have an identifier, it's completed by the whitespace
                    if !identifier.is_empty() { break 'ident; }
                } else {
                    // Consume c
                    consumed_amount += c.len_utf8();
                    //  and add it to the identifier
                    identifier.push(c);
                }
            }
        }
        // Consume character(s) from final loop
        if consumed_amount != 0 {
            read.consume(consumed_amount);
        }
        // Identifiers cannot be empty
        if identifier.is_empty() {
            return Ok(None);
        }
        // :tada:
        return Ok(Some(identifier));
    }
}

impl fmt::Display for Graph {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "<!-- graph nodes -->")?;
        for node in &self.nodes {
            writeln!(f, r#"<text font-style="italic" text-anchor="end" x="{}" y="{}">{}</text>"#,
                node.position[0] - 8., node.position[1], node.name)?;
            if node.draw_circle {
                writeln!(f, r#"<circle cx="{}" cy="{}" r="2" />"#,
                    node.position[0], node.position[1] - 5.)?;
            }
        }
        writeln!(f, "<!-- graph edges -->")?;
        for (from, destinations) in self.edges.iter().enumerate() {
            for &to in destinations {
                let from = &self.nodes[from];
                let to = &self.nodes[to];
                writeln!(f, r#"<line x1="{}" y1="{}" x2="{}" y2="{}" marker-end="url(#triangle)" stroke="black" />"#,
                    from.position[0], from.position[1] - 5.,
                    to.position[0],   to.position[1]   - 5.)?;
            }
        }
        Ok(())
    }
}

pub enum InputError {
    IoError(io::Error),
    /// catch all error for bad inputs
    ParseError,
    /// unexpected EOF
    ShortRead,
    /// couldn't read a line as a number
    CantParseNumber(String, ParseIntError),
    /// couldn't read a line as a floating point number
    CantParseFloat(String, ParseFloatError),
    /// each node needs a unique name
    NodesHaveSameName(String),
    /// edge has invalid node
    EdgeReferencesNonexistantNode(String),
}

impl fmt::Debug for InputError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        use InputError::*;
        match self {
            ParseError => write!(f, "error parsing"),
            ShortRead => write!(f, "error parsing: short read"),
            IoError(e) => write!(f, "error reading: {}", e),
            CantParseNumber(s, e) => write!(f, "error parsing {:?} as an integer: {}", s, e),
            CantParseFloat(s, e) => write!(f, "error parsing {:?} as a decimal number: {}", s, e),
            NodesHaveSameName(s) => write!(f, "two nodes have name {:?}", s),
            EdgeReferencesNonexistantNode(s) => write!(f, "edge references nonexistant node {:?}", s),
        }
    }
}

impl From<io::Error> for InputError {
    fn from(e: io::Error) -> Self {
        InputError::IoError(e)
    }
}

/// determine the level of every node
///
/// (Implements ROUTINE LABEL.1)
fn determine_node_level(graph: &Graph) -> (Vec<Vec<usize>>, Vec<usize>) {
    let mut level = vec![usize::MAX; graph.nodes.len()];
    let mut nodes_in_level = vec![];
    let mut explore_stack = vec![];
    for node in 0..graph.nodes.len() {
        // skip nodes which have a level assigned
        if level[node] != usize::MAX {
            continue;
        }
        
        // visit node and its descendants (with indeterminate level) in a
        // pre-order traversal, pushing them onto the reverse_queue for a
        // post-order traversal.
        explore_stack.push(node);
        while let Some(node) = explore_stack.pop() {
            // max is the level of this node
            // before we look at our children, assume we're in level 0.
            let mut max = 0;
            let mut exploring = false;
            for &child in &graph.edges[node] {
                if exploring {
                    if level[child] == usize::MAX {
                        explore_stack.push(child);
                    }
                } else {
                    if level[child] == usize::MAX {
                        // we haven't visited a child of this node yet.
                        // visit all unvisited descendants of this node and then revisit this node.
                        exploring = true;
                        explore_stack.push(node);
                        explore_stack.push(child);
                    } else if level[child] + 1 > max {
                        // we must be in a larger level than our child.
                        max = level[child] + 1;
                    }
                }
            }
            if !exploring {
                // ensure nodes_in_level has a vector for our level.
                if nodes_in_level.len() <= max {
                    nodes_in_level.resize(max+1, vec![]);
                }
                // insert ourselves in both mappings
                nodes_in_level[max].push(node);
                level[node] = max;
            }
        }
    }

    (nodes_in_level, level)
}

fn main() {
    let mut text = br#"
    MLDAG
    // node count
    17
    // name, posx, posy
    A 100 0
    B 200 0
    C 300 0
    D 180 0
    E 250 0
    F 350 0
    G 80  0
    H 200 0
    I 300 0
    J 400 0
    K 40  0
    L 160 0
    M 150 0
    N 200 0
    O 100 0
    P 200 0
    Q 300 0
    // edge count
    23
    // from, to
    O N
    P N
    Q N
    Q J
    N M
    N I
    M K
    M L
    K G
    L G
    L H
    G A
    G E
    H D
    H F
    I D
    I F
    J F
    D B
    D C
    E C
    F A
    F C
    "# as &[u8];

    let graph = Graph::try_parse(&mut text);
    let mut graph = match graph {
        Ok(graph) => graph,
        Err(e) => {
            eprintln!("Couldn't parse input DAG: {:?}", e);
            return;
        }
    };

    
    // set their Y to their level
    let (nodes_in_level, level) = determine_node_level(&graph);
    for (nidx, node) in graph.nodes.iter_mut().enumerate() {
        node.position[1] = (nodes_in_level.len() - level[nidx]) as f32 * LEVEL_HEIGHT;
        node.draw_circle = nodes_in_level.len() != level[nidx]+1;
    }

    println!("<!DOCTYPE html><html><title>Lexicographical Ordering in Linear Time</title>");
    let label = routine_label(&graph);
    
    svg_start("FINAL RESULT");
    println!("<!-- routine labels -->");
    //for (node, label) in level.iter().enumerate() {
    for (node, label) in label.iter().enumerate() {
        let node = &graph.nodes[node];
        println!(r#"<text font-style="italic" text-anchor="start" x="{}" y="{}">{}</text>"#,
            node.position[0] + 8., node.position[1], label+1);
    }
    svg_end(&graph);

    println!("{}", r#"
    <style>
    svg:not(:last-of-type) { display:block; margin: 0 auto 100px auto; }
    svg:last-of-type { height: 100vh; display: block; margin: 0 auto; }
    html, body {margin: 0; padding: 0;}
    </style>
    <script>
    var active_slide = document.getElementsByTagName("svg")[0];
    active_slide.scrollIntoView();
    document.addEventListener("keydown", function(e){
        var next;
        if (e.code == "Home") {
            next = document.getElementsByTagName("svg")[0];
        } else if (e.code == "End") {
            next = document.getElementsByTagName("svg");
            next = next[next.length - 1]
        } else if (e.code == "ArrowRight" || e.code == "KeyD") {
            next = active_slide.nextElementSibling;
        } else if (e.code == "ArrowLeft" || e.code == "KeyA") {
            next = active_slide.previousElementSibling;
        } else {
            return;
        }
        if (!next || next.tagName != "svg") {
            return;
        }
        active_slide = next;
        active_slide.scrollIntoView();
        return true;
    })
    </script>
    </html>"#);
}

fn svg_start( comment: &str ){
    println!(r#"<svg width="{}" height="400">
    <defs>
        <marker id="triangle" viewBox="0 0 10 10"
          refX="12" refY="5"
          markerUnits="strokeWidth"
          markerWidth="10" markerHeight="10"
          orient="auto">
            <path d="M 4 2 L 10 5 L 4 8 z"/>
        </marker>
        <pattern id="hazard" width="25" height="25" patternUnits="userSpaceOnUse">
            <rect x="0" y="0" width="25" height="25" fill="black" opacity="0.5"/>
            <line x1="0" y1="-12.5" x2="37.5" y2="25" stroke-width="10" stroke="black" />
            <line x1="-12.5" y1="0" x2="25" y2="37.5" stroke-width="10" stroke="black" />
        </pattern>
    </defs>
    <text text-anchor="middle" x="{}" y="20" font-weight="bold">{}</text>
    "#, GRAPH_WIDTH, GRAPH_WIDTH/2., comment);
}

fn svg_end( graph: &Graph ){
    println!("{}</svg>", graph);
}
