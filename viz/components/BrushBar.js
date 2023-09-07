import React, { useEffect, useRef } from 'react';
import * as d3 from 'd3';

const BrushingBarChart = ({ data , width, height , handleFilter}) => {
  const svgRef = useRef(null);

  useEffect(() => {

    const svg = d3.select(svgRef.current);

    // Define margins and inner dimensions
    const margin = { top: 20, right: 20, bottom: 20, left: 20 };
    const innerWidth = width - margin.left - margin.right;
    const innerHeight = height - margin.top - margin.bottom;

    // Create scale for the x-axis
    const xScale = d3.scaleLinear()
      .domain([d3.min(data), d3.max(data)])
      .range([0, innerWidth]);

    // Create histogram bins
    const histogram = d3.histogram()
      .domain(xScale.domain())
      .thresholds(xScale.ticks(10));

    const bins = histogram(data);

    // Create scale for the y-axis
    const yScale = d3.scaleLinear()
      .domain([0, d3.max(bins, d => d.length)])
      .range([innerHeight, 0]);

    // Create brush
    const brush = d3.brushX()
      .extent([[0, 0], [innerWidth, innerHeight]])
      .on('brush', function(event){
        // brushGroup.remove()
        // brushGroup.call(brush.move, null);
        d3.select(".brush").call(brush.extent([0, 0]))
        if (event.selection) {
          const [x0, x1] = event.selection;
          const selectedBins = bins.filter(bin =>
            xScale(bin.x0) >= x0 && xScale(bin.x1) <= x1
          );

          svg.selectAll('.bar')
            .style('fill', d =>
              selectedBins.includes(d) ? 'steelblue' : 'gray'
            );
            // d3.select(this).call(brushHandle, event.selection);
        }else{
          svg.selectAll('.bar')
            .style('fill', 'steelblue'
            );
        }
      })
      .on('end', (event) => {
        if (event.selection) {
          //update bar color 
          const [x0, x1] = event.selection;
          const selectedBins = bins.filter(bin =>
            xScale(bin.x0) >= x0 && xScale(bin.x1) <= x1
          );
            console.log(selectedBins)
          svg.selectAll('.bar')
            .style('fill', d =>
              selectedBins.includes(d) ? 'steelblue' : 'gray'
            );
            //filter dataset 
            const minbin = Math.min(...[].concat(...selectedBins))
            const maxbin = Math.max(...[].concat(...selectedBins))
            handleFilter([minbin,maxbin])
        }else{
          svg.selectAll('.bar')
          .style('fill', 'steelblue'
          );

        }
      });

    // Create chart container
    const chart = svg.append('g')
      .attr('transform', `translate(${margin.left},${margin.top})`);

    // Draw histogram bars
    chart.selectAll('.bar')
      .data(bins)
      .enter()
      .append('rect')
      .attr('class', 'bar')
      .attr('x', d => xScale(d.x0))
      .attr('y', d => yScale(d.length))
      .attr('width', d => xScale(d.x1) - xScale(d.x0) - 1)
      .attr('height', d => innerHeight - yScale(d.length))
      .style('fill', 'gray');

    // Add x-axis
    const xAxis = d3.axisBottom(xScale);
    chart.append('g')
      .attr('transform', `translate(0,${innerHeight})`)
      .call(xAxis);

    // Add y-axis
    const yAxisTicks = yScale.ticks()
    .filter(tick => Number.isInteger(tick));

    const yAxis = d3.axisLeft(yScale).tickValues(yAxisTicks).tickFormat(d3.format('d'));
    chart.append('g')
      .call(yAxis);

    // Add brush to the chart
    const brushGroup = chart.append('g')
      .attr('class', 'brush')
      .call(brush);

    // Add handles to the brush
    // const brushHandle = (g, s) => g.selectAll('.handle--custom')
    //   .data([{ type: 'w' }, { type: 'e' }])
    //   .join(
    //     enter => enter.append("path")
    //         // .attr("class", "handle--custom")
    //         .attr('class', 'handle--custom')
    //         .attr('cursor', 'ew-resize')
    //         .attr("fill", "#eee")
    //         .attr("fill-opacity", 0.8)
    //         .attr("stroke", "#666")
    //         .attr("stroke-width", 1.5)
    //         .attr('d', d => {
    //           const e = +(d.type === 'e');
    //           const x = e ? 1 : -1;
    //           const y = innerHeight / 2;
      
    //           return `M${0.5 * x},${y}A6,6 0 0 ${e} ${6.5 * x},${y + 6}V${2 * y - 6}A6,6 0 0 ${e} ${0.5 * x},${2 * y}ZM${2.5 * x},${y + 8}V${2 * y - 8}M${4.5 * x},${y + 8}V${2 * y - 8}`;
    //         })
    //   )
    //   .attr('transform', d => `translate(${xScale.range()[d.type === 'e' ? 1 : 0]},0)`);


  }, [data, width, height]);

  return (
    <svg ref={svgRef} width={width} height={height}></svg>
  );
};

export default BrushingBarChart;