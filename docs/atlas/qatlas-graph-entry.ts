// qatlas-graph-entry.ts — a standalone, data-driven port of Quartz's
// graph.inline.ts (d3-force + PixiJS) for embedding in Documenter.jl (or any
// site).  Stripped of Quartz internals (contentIndex/slugs/SPA-nav/tags/depth):
// it takes {nodes, links} directly and renders the Obsidian-style force graph.
//
// Encoding:
//   node color  = node type (model / quantity / universality / bound)
//   edge color  = verified status of the result (exact & universal share a
//                 color; bound; approx); `realizes` edges use a neutral color
//   edge style  = solid (verified) / dashed (not verified)
//   gaps are NOT marked explicitly — an isolated node *is* the gap.
//
// Bundle to an IIFE that sets window.QAtlasGraph:
//   esbuild qatlas-graph-entry.ts --bundle --format=iife --minify \
//     --outfile=qatlas-graph.js
//   window.QAtlasGraph.render(container, {nodes, links}, cfg)

import {
  SimulationNodeDatum,
  SimulationLinkDatum,
  Simulation,
  forceSimulation,
  forceManyBody,
  forceCenter,
  forceLink,
  forceCollide,
  zoomIdentity,
  select,
  drag,
  zoom,
} from "d3"
import { Text, Graphics, Application, Container, Circle } from "pixi.js"
import { Group as TweenGroup, Tween as Tweened } from "@tweenjs/tween.js"

type NodeIn = { id: string; text?: string; group?: string; url?: string }
type LinkIn = {
  source: string
  target: string
  kind?: string // "realizes" | "provides"
  status?: string // exact | bound | approx | universal
  verified?: boolean
}

type NodeData = {
  id: string
  text: string
  group: string
  url?: string
} & SimulationNodeDatum

type LinkData = {
  source: NodeData
  target: NodeData
  kind: string
  status: string
  verified: boolean
} & SimulationLinkDatum<NodeData>

type GraphicsInfo = { color: number; gfx: Graphics; alpha: number; active: boolean }
type LinkRenderData = GraphicsInfo & { simulationData: LinkData; verified: boolean }
type NodeRenderData = GraphicsInfo & { simulationData: NodeData; label: Text }

type Cfg = {
  drag?: boolean
  zoom?: boolean
  scale?: number
  repelForce?: number
  centerForce?: number
  linkDistance?: number
  linkStrength?: number
  fontSize?: number
  opacityScale?: number
  focusOnHover?: boolean
  legend?: boolean
  search?: boolean
  colors?: Record<string, string>
  statusColors?: Record<string, string>
  realizesColor?: string
  realizesLabel?: string
  labelColor?: string
  labelStroke?: string
}

const DEFAULT_COLORS: Record<string, string> = {
  // Model is red/pink so its (many) nodes don't blend into the blue
  // exact/universal edges; Bound domain takes the freed-up blue.
  model: "#c95f5f",
  class: "#9b6cc9",
  bound: "#4f8cc9",
  quantity: "#5fae6f",
  default: "#888888",
}

// exact & universal deliberately share a color (low visual load).
const STATUS_COLORS: Record<string, string> = {
  exact: "#6f9fd8",
  universal: "#6f9fd8",
  bound: "#d88f8f",
  approx: "#e0a84d",
  default: "#9aa0a8",
}

function hexToNum(hex: string): number {
  return parseInt(hex.replace("#", ""), 16)
}

type TweenNode = { update: (t: number) => void; stop: () => void }

async function renderGraph(
  graph: HTMLElement,
  dataIn: { nodes: NodeIn[]; links: LinkIn[] },
  cfg: Cfg = {},
) {
  const {
    drag: enableDrag = true,
    zoom: enableZoom = true,
    scale = 1.0,
    repelForce = 0.85,
    centerForce = 0.2,
    linkDistance = 105,
    linkStrength = 0.07,
    fontSize = 0.6,
    opacityScale = 1.1,
    focusOnHover = true,
    colors = {},
    statusColors = {},
    realizesColor = "#9aa0a8",
    realizesLabel = "realizes",
    labelColor = "#eaeaea",
    labelStroke = "#15171a",
  } = cfg
  const palette = { ...DEFAULT_COLORS, ...colors }
  const statusPalette = { ...STATUS_COLORS, ...statusColors }
  const colorOf = (g: string) => hexToNum(palette[g] ?? palette.default)
  const linkColorOf = (l: LinkData) =>
    hexToNum(
      l.kind === "realizes" ? realizesColor : (statusPalette[l.status] ?? statusPalette.default),
    )

  while (graph.firstChild) graph.removeChild(graph.firstChild)

  const nodes: NodeData[] = dataIn.nodes.map((n) => ({
    id: n.id,
    text: n.text ?? n.id,
    group: n.group ?? "default",
    url: n.url,
  }))
  const byId = new Map(nodes.map((n) => [n.id, n]))
  const links: LinkData[] = dataIn.links
    .filter((l) => byId.has(l.source) && byId.has(l.target))
    .map((l) => ({
      source: byId.get(l.source)!,
      target: byId.get(l.target)!,
      kind: l.kind ?? "provides",
      status: l.status ?? "exact",
      verified: l.verified ?? false,
    }))
  const graphData = { nodes, links }

  const width = graph.offsetWidth || 800
  const height = Math.max(graph.offsetHeight, 250)

  const simulation: Simulation<NodeData, LinkData> = forceSimulation<NodeData>(graphData.nodes)
    .force("charge", forceManyBody().strength(-100 * repelForce))
    .force("center", forceCenter().strength(centerForce))
    .force("link", forceLink(graphData.links).distance(linkDistance).strength(linkStrength))
    .force("collide", forceCollide<NodeData>((n) => nodeRadius(n) + 6).iterations(3))

  function nodeRadius(d: NodeData) {
    const numLinks = graphData.links.filter(
      (l) => l.source.id === d.id || l.target.id === d.id,
    ).length
    return 3 + Math.sqrt(numLinks)
  }

  let hoveredNodeId: string | null = null
  let hoveredNeighbours: Set<string> = new Set()
  let searchMatches: Set<string> = new Set()
  const linkRenderData: LinkRenderData[] = []
  const nodeRenderData: NodeRenderData[] = []
  const tweens = new Map<string, TweenNode>()

  function updateHoverInfo(newId: string | null) {
    hoveredNodeId = newId
    if (newId === null) {
      hoveredNeighbours = new Set()
      for (const n of nodeRenderData) n.active = false
      for (const l of linkRenderData) l.active = false
    } else {
      hoveredNeighbours = new Set()
      for (const l of linkRenderData) {
        const d = l.simulationData
        if (d.source.id === newId || d.target.id === newId) {
          hoveredNeighbours.add(d.source.id)
          hoveredNeighbours.add(d.target.id)
        }
        l.active = d.source.id === newId || d.target.id === newId
      }
      for (const n of nodeRenderData) n.active = hoveredNeighbours.has(n.simulationData.id)
    }
  }

  let dragStartTime = 0
  let dragging = false

  function renderLinks() {
    tweens.get("link")?.stop()
    const tg = new TweenGroup()
    for (const l of linkRenderData) {
      let alpha = 1
      if (hoveredNodeId) alpha = l.active ? 1 : 0.18
      else if (searchMatches.size) {
        const d = l.simulationData
        const s = searchMatches.has(d.source.id),
          t = searchMatches.has(d.target.id)
        alpha = s && t ? 1 : s || t ? 0.4 : 0.05
      }
      l.color = linkColorOf(l.simulationData)
      tg.add(new Tweened<LinkRenderData>(l).to({ alpha }, 200))
    }
    tg.getAll().forEach((t) => t.start())
    tweens.set("link", {
      update: tg.update.bind(tg),
      stop: () => tg.getAll().forEach((t) => t.stop()),
    })
  }

  function renderLabels() {
    tweens.get("label")?.stop()
    const tg = new TweenGroup()
    const def = 1 / scale
    const act = def * 1.1
    const focusing = hoveredNodeId !== null || searchMatches.size > 0
    for (const n of nodeRenderData) {
      const id = n.simulationData.id
      const isH = hoveredNodeId === id
      const show = isH || hoveredNeighbours.has(id) || searchMatches.has(id)
      const alpha = show ? 1 : focusing ? 0 : n.label.alpha
      tg.add(
        new Tweened<Text>(n.label).to(
          { alpha, scale: { x: isH ? act : def, y: isH ? act : def } },
          100,
        ),
      )
    }
    tg.getAll().forEach((t) => t.start())
    tweens.set("label", {
      update: tg.update.bind(tg),
      stop: () => tg.getAll().forEach((t) => t.stop()),
    })
  }

  function renderNodes() {
    tweens.get("hover")?.stop()
    const tg = new TweenGroup()
    for (const n of nodeRenderData) {
      let alpha = 1
      if (hoveredNodeId !== null && focusOnHover) alpha = n.active ? 1 : 0.2
      else if (searchMatches.size) alpha = searchMatches.has(n.simulationData.id) ? 1 : 0.12
      tg.add(new Tweened<Graphics>(n.gfx, tg).to({ alpha }, 200))
    }
    tg.getAll().forEach((t) => t.start())
    tweens.set("hover", {
      update: tg.update.bind(tg),
      stop: () => tg.getAll().forEach((t) => t.stop()),
    })
  }

  function renderPixiFromD3() {
    renderNodes()
    renderLinks()
    renderLabels()
  }

  const app = new Application()
  await app.init({
    width,
    height,
    antialias: true,
    autoStart: false,
    autoDensity: true,
    backgroundAlpha: 0,
    resolution: window.devicePixelRatio,
    eventMode: "static",
  })
  graph.appendChild(app.canvas)

  if (cfg.legend !== false) {
    graph.style.position = graph.style.position || "relative"
    const present = new Set(nodes.map((n) => n.group))
    const presentStatus = new Set(links.filter((l) => l.kind !== "realizes").map((l) => l.status))
    const hasRealizes = links.some((l) => l.kind === "realizes")
    const nodeItems: [string, string][] = [
      ["model", "Model"],
      ["class", "Universality class"],
      ["bound", "Bound domain"],
      ["quantity", "Quantity"],
    ]
    // exact & universal collapse to one swatch
    const statusItems: [string, string, string][] = [
      ["exact", STATUS_COLORS.exact, "exact / universal"],
      ["bound", STATUS_COLORS.bound, "bound"],
      ["approx", STATUS_COLORS.approx, "approx"],
    ]
    const dot = (c: string) =>
      `<span style="display:inline-block;width:11px;height:11px;border-radius:50%;background:${c};margin-right:7px;vertical-align:middle"></span>`
    const solid = (c: string) =>
      `<span style="display:inline-block;width:17px;border-top:2px solid ${c};margin:0 7px 4px 0;vertical-align:middle"></span>`
    const dashed = (c: string) =>
      `<span style="display:inline-block;width:17px;border-top:2px dashed ${c};margin:0 7px 4px 0;vertical-align:middle"></span>`
    const nodeHtml = nodeItems
      .filter(([g]) => present.has(g))
      .map(([g, label]) => `<div>${dot(palette[g])}${label}</div>`)
      .join("")
    const statusHtml = statusItems
      .filter(([k]) => presentStatus.has(k) || (k === "exact" && presentStatus.has("universal")))
      .map(([, c, label]) => `<div>${solid(c)}${label}</div>`)
      .join("")
    const realizesHtml = hasRealizes ? `<div>${solid(realizesColor)}${realizesLabel}</div>` : ""
    const styleHtml =
      `<div style="margin-top:3px">${solid("#cfcfcf")}verified</div>` +
      `<div>${dashed("#cfcfcf")}not verified</div>`
    const sep = (t: string) => `<div style="margin:5px 0 2px;opacity:0.5">${t}</div>`
    const legend = document.createElement("div")
    legend.style.cssText =
      "position:absolute;top:10px;left:10px;padding:8px 11px;font-size:12px;line-height:1.7;" +
      "background:rgba(18,20,24,0.6);color:#e8e8e8;border-radius:6px;pointer-events:none"
    legend.innerHTML =
      nodeHtml + sep("edges") + statusHtml + realizesHtml + sep("style") + styleHtml
    graph.appendChild(legend)
  }

  // Search box: filter nodes by name, highlighting matches (and their labels)
  // while dimming the rest — the graph has no other search affordance.
  if (cfg.search !== false) {
    graph.style.position = graph.style.position || "relative"
    const box = document.createElement("input")
    box.type = "search"
    box.placeholder = "search nodes…"
    box.style.cssText =
      "position:absolute;top:10px;right:10px;z-index:5;width:150px;padding:5px 9px;" +
      "font-size:12px;background:rgba(18,20,24,0.85);color:#e8e8e8;" +
      "border:1px solid #3a3f47;border-radius:6px;outline:none"
    box.addEventListener("input", () => {
      const q = box.value.trim().toLowerCase()
      searchMatches = new Set()
      if (q)
        for (const n of nodes)
          if (n.text.toLowerCase().includes(q) || n.id.toLowerCase().includes(q))
            searchMatches.add(n.id)
      renderPixiFromD3()
    })
    graph.appendChild(box)
  }

  const stage = app.stage
  stage.interactive = false
  const labelsC = new Container<Text>({ zIndex: 3, isRenderGroup: true })
  const nodesC = new Container<Graphics>({ zIndex: 2, isRenderGroup: true })
  const linkC = new Container<Graphics>({ zIndex: 1, isRenderGroup: true })
  stage.addChild(nodesC, labelsC, linkC)

  for (const n of graphData.nodes) {
    const label = new Text({
      interactive: false,
      eventMode: "none",
      text: n.text,
      alpha: 0,
      anchor: { x: 0.5, y: 1.2 },
      style: {
        fontSize: fontSize * 16,
        fill: hexToNum(labelColor),
        stroke: { color: hexToNum(labelStroke), width: 3 },
      },
      resolution: window.devicePixelRatio * 2,
    })
    label.scale.set(1 / scale)
    let oldOpacity = 0
    const gfx = new Graphics({
      interactive: true,
      label: n.id,
      eventMode: "static",
      hitArea: new Circle(0, 0, nodeRadius(n)),
      cursor: "pointer",
    })
      .circle(0, 0, nodeRadius(n))
      .fill({ color: colorOf(n.group) })
      .on("pointerover", (e: any) => {
        updateHoverInfo(e.target.label)
        oldOpacity = label.alpha
        if (!dragging) renderPixiFromD3()
      })
      .on("pointerleave", () => {
        updateHoverInfo(null)
        label.alpha = oldOpacity
        if (!dragging) renderPixiFromD3()
      })
    nodesC.addChild(gfx)
    labelsC.addChild(label)
    nodeRenderData.push({
      simulationData: n,
      gfx,
      label,
      color: colorOf(n.group),
      alpha: 1,
      active: false,
    })
  }

  for (const l of graphData.links) {
    const gfx = new Graphics({ interactive: false, eventMode: "none" })
    linkC.addChild(gfx)
    linkRenderData.push({
      simulationData: l,
      gfx,
      verified: l.verified,
      color: linkColorOf(l),
      alpha: 1,
      active: false,
    })
  }

  let currentTransform = zoomIdentity
  function navTo(n: NodeData) {
    if (n.url) window.location.href = n.url
  }
  if (enableDrag) {
    select<HTMLCanvasElement, NodeData | undefined>(app.canvas as any).call(
      drag<HTMLCanvasElement, NodeData | undefined>()
        .container(() => app.canvas as any)
        .subject(() => graphData.nodes.find((n) => n.id === hoveredNodeId))
        .on("start", (event: any) => {
          if (!event.active) simulation.alphaTarget(1).restart()
          event.subject.fx = event.subject.x
          event.subject.fy = event.subject.y
          event.subject.__p = { x: event.subject.x, y: event.subject.y }
          dragStartTime = Date.now()
          dragging = true
        })
        .on("drag", (event: any) => {
          const p = event.subject.__p
          event.subject.fx = p.x + (event.x - p.x) / currentTransform.k
          event.subject.fy = p.y + (event.y - p.y) / currentTransform.k
        })
        .on("end", (event: any) => {
          if (!event.active) simulation.alphaTarget(0)
          event.subject.fx = null
          event.subject.fy = null
          dragging = false
          if (Date.now() - dragStartTime < 500) {
            const node = graphData.nodes.find((n) => n.id === event.subject.id)
            if (node) navTo(node)
          }
        }),
    )
  } else {
    for (const nd of nodeRenderData) nd.gfx.on("click", () => navTo(nd.simulationData))
  }

  if (enableZoom) {
    select<HTMLCanvasElement, NodeData>(app.canvas as any).call(
      zoom<HTMLCanvasElement, NodeData>()
        .extent([
          [0, 0],
          [width, height],
        ])
        .scaleExtent([0.25, 4])
        .on("zoom", ({ transform }: any) => {
          currentTransform = transform
          stage.scale.set(transform.k, transform.k)
          stage.position.set(transform.x, transform.y)
          const s = transform.k * opacityScale
          const sOp = Math.max((s - 1) / 3.75, 0)
          const actLabels = nodeRenderData.filter((n) => n.active).map((n) => n.label)
          for (const lab of labelsC.children) if (!actLabels.includes(lab)) lab.alpha = sOp
        }),
    )
  }

  let stop = false
  function animate(time: number) {
    if (stop) return
    for (const n of nodeRenderData) {
      const { x, y } = n.simulationData
      if (x == null || y == null) continue
      n.gfx.position.set(x + width / 2, y + height / 2)
      n.label.position.set(x + width / 2, y + height / 2)
    }
    for (const l of linkRenderData) {
      const d = l.simulationData
      const sx = d.source.x! + width / 2,
        sy = d.source.y! + height / 2
      const tx = d.target.x! + width / 2,
        ty = d.target.y! + height / 2
      const dx = tx - sx,
        dy = ty - sy
      const len = Math.hypot(dx, dy) || 1
      const ux = dx / len,
        uy = dy / len
      const tr = nodeRadius(d.target) + 1.5
      const end = Math.max(len - tr, 0) // stop at the target node's edge
      const ex = sx + ux * end,
        ey = sy + uy * end
      l.gfx.clear()
      if (l.verified) {
        l.gfx.moveTo(sx, sy).lineTo(ex, ey).stroke({ alpha: l.alpha, width: 1, color: l.color })
      } else {
        // dashed line for unverified edges
        const dash = 5,
          gp = 4
        let dd = 0
        while (dd < end) {
          const x1 = sx + ux * dd,
            y1 = sy + uy * dd
          const d2 = Math.min(dd + dash, end)
          l.gfx.moveTo(x1, y1).lineTo(sx + ux * d2, sy + uy * d2)
          dd += dash + gp
        }
        l.gfx.stroke({ alpha: l.alpha, width: 1, color: l.color })
      }
      // arrowhead at the target end
      const ah = 4.5,
        a = Math.atan2(dy, dx)
      l.gfx
        .moveTo(ex, ey)
        .lineTo(ex - ah * Math.cos(a - 0.45), ey - ah * Math.sin(a - 0.45))
        .lineTo(ex - ah * Math.cos(a + 0.45), ey - ah * Math.sin(a + 0.45))
        .lineTo(ex, ey)
        .fill({ alpha: l.alpha, color: l.color })
    }
    tweens.forEach((t) => t.update(time))
    app.renderer.render(stage)
    requestAnimationFrame(animate)
  }
  requestAnimationFrame(animate)
  return () => {
    stop = true
    app.destroy()
  }
}

;(window as any).QAtlasGraph = {
  render(el: HTMLElement | string, data: { nodes: NodeIn[]; links: LinkIn[] }, cfg: Cfg = {}) {
    const node = typeof el === "string" ? document.getElementById(el) : el
    if (!node) return
    renderGraph(node as HTMLElement, data, cfg).catch((err) => {
      ;(node as HTMLElement).innerHTML =
        "<p style='padding:1em'>graph render error: " +
        (err && err.message ? err.message : err) +
        "</p>"
    })
  },
}
