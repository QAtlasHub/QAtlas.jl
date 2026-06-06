// qatlas-graph-entry.ts — a standalone, data-driven port of Quartz's
// graph.inline.ts (d3-force + PixiJS) for embedding in Documenter.jl (or any
// site).  Stripped of Quartz internals (contentIndex/slugs/SPA-nav/tags/depth):
// it takes {nodes, links} directly and renders the Obsidian-style force graph.
//
// Bundle to an IIFE that sets window.QAtlasGraph:
//   esbuild qatlas-graph-entry.ts --bundle --format=iife --global-name=QAtlasGraph \
//     --minify --outfile=qatlas-graph.js
//
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
type LinkIn = { source: string; target: string; gap?: boolean; label?: string }

type NodeData = {
  id: string
  text: string
  group: string
  url?: string
} & SimulationNodeDatum

type LinkData = {
  source: NodeData
  target: NodeData
  gap: boolean
} & SimulationLinkDatum<NodeData>

type GraphicsInfo = { color: number; gfx: Graphics; alpha: number; active: boolean }
type LinkRenderData = GraphicsInfo & { simulationData: LinkData; gap: boolean }
type NodeRenderData = GraphicsInfo & { simulationData: NodeData; label: Text }

type Cfg = {
  drag?: boolean
  zoom?: boolean
  scale?: number
  repelForce?: number
  centerForce?: number
  linkDistance?: number
  fontSize?: number
  opacityScale?: number
  focusOnHover?: boolean
  legend?: boolean
  colors?: Record<string, string>
  linkColor?: string
  gapColor?: string
  labelColor?: string
}

const DEFAULT_COLORS: Record<string, string> = {
  model: "#4f8cc9",
  class: "#9b6cc9",
  bound: "#c95f5f",
  quantity: "#5fae6f",
  gap: "#e8a33d",
  default: "#888888",
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
    repelForce = 0.5,
    centerForce = 0.3,
    linkDistance = 30,
    fontSize = 0.5,
    opacityScale = 1.0,
    focusOnHover = true,
    colors = {},
    linkColor = "#b8b8b8",
    gapColor = "#e8a33d",
    labelColor = "#2b2b2b",
  } = cfg
  const palette = { ...DEFAULT_COLORS, ...colors }
  const colorOf = (g: string) => hexToNum(palette[g] ?? palette.default)

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
    .map((l) => ({ source: byId.get(l.source)!, target: byId.get(l.target)!, gap: !!l.gap }))
  const graphData = { nodes, links }

  const width = graph.offsetWidth || 800
  const height = Math.max(graph.offsetHeight, 250)

  const simulation: Simulation<NodeData, LinkData> = forceSimulation<NodeData>(graphData.nodes)
    .force("charge", forceManyBody().strength(-100 * repelForce))
    .force("center", forceCenter().strength(centerForce))
    .force("link", forceLink(graphData.links).distance(linkDistance))
    .force("collide", forceCollide<NodeData>((n) => nodeRadius(n)).iterations(3))

  function nodeRadius(d: NodeData) {
    const numLinks = graphData.links.filter(
      (l) => l.source.id === d.id || l.target.id === d.id,
    ).length
    return 3 + Math.sqrt(numLinks)
  }

  let hoveredNodeId: string | null = null
  let hoveredNeighbours: Set<string> = new Set()
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
      if (hoveredNodeId) alpha = l.active ? 1 : 0.2
      l.color = l.gap ? hexToNum(gapColor) : hexToNum(linkColor)
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
    for (const n of nodeRenderData) {
      const id = n.simulationData.id
      const isH = hoveredNodeId === id
      const show = isH || hoveredNeighbours.has(id)
      const alpha = show ? 1 : hoveredNodeId !== null ? 0 : n.label.alpha
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
    const legendItems: [string, string][] = [
      ["model", "Model"],
      ["class", "Universality class"],
      ["bound", "Bound domain"],
      ["quantity", "Quantity"],
      ["gap", "Coherence gap"],
    ]
    const legend = document.createElement("div")
    legend.style.cssText =
      "position:absolute;top:10px;left:10px;padding:7px 10px;font-size:12px;line-height:1.7;" +
      "background:rgba(127,127,127,0.14);border-radius:6px;pointer-events:none"
    legend.innerHTML = legendItems
      .filter(([g]) => present.has(g))
      .map(
        ([g, label]) =>
          `<div><span style="display:inline-block;width:11px;height:11px;border-radius:50%;` +
          `background:${palette[g]};margin-right:7px;vertical-align:middle"></span>${label}</div>`,
      )
      .join("")
    graph.appendChild(legend)
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
      style: { fontSize: fontSize * 16, fill: hexToNum(labelColor) },
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
      gap: l.gap,
      color: l.gap ? hexToNum(gapColor) : hexToNum(linkColor),
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
      const ex = tx - ux * tr,
        ey = ty - uy * tr // line stops at the target node's edge
      l.gfx.clear()
      l.gfx.moveTo(sx, sy).lineTo(ex, ey).stroke({ alpha: l.alpha, width: 1, color: l.color })
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
