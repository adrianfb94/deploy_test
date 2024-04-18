const Logger = require("discord-blackhole.logger");
const mkdirp = require("mkdirp-sync");
const fs = require("fs");

class myLogger extends Logger.class
{
	constructor (path)
	{
		super();
		this.path = path;
		this.assureFile(this.path);
	}
	
	assureDirectory (dirname)
	{
		if (!fs.existsSync(dirname)) mkdirp(dirname);
	}
	
	assureFile (filename)
	{
		if (!fs.existsSync(require("path").dirname(filename))) this.assureDirectory(require("path").dirname(filename));
		if (!fs.existsSync(filename)) fs.writeFileSync(filename, "");
	}
	
	append (content)
	{
		this.assureFile(this.path);
		fs.appendFile(this.path, content, function (err) {
			if (err) process.stdout.write("ERROR Failed to write file!");
		});
	}
	
	timeNoColor ()
	{
		var time = this.format;
		const Time = this.getTimeObject();
		for (var key in Time)
		{
			const value = Time[key];
			time = time.split("{" + key + "}").join(value);
		}
		
		return time;
	}
	
	formatNoColor (type, message)
	{
		const v = (this.timeNoColor() + " " + type + " " + message.join(" ") + "\n").replace(/\u001B\[[0-9][0-9]m/gi, "");
		return v;
	}
	
	log ()
	{
		const args = Array.prototype.slice.call(arguments);
		const msg = this.formatNoColor("   ", args);
		this.append(msg);
		super.log(...args);
		return this;
	}
	
	debug ()
	{
		const args = Array.prototype.slice.call(arguments);
		const msg = this.formatNoColor("DEB", args);
		this.append(msg);
		super.debug(...args);
		return this;
	}
	
	info ()
	{
		const args = Array.prototype.slice.call(arguments);
		const msg = this.formatNoColor("INF", args);
		this.append(msg);
		super.info(...args);
		return this;
	}
	
	alert ()
	{
		const args = Array.prototype.slice.call(arguments);
		const msg = this.formatNoColor("ALE", args);
		this.append(msg);
		super.alert(...args);
		return this;
	}
	
	success ()
	{
		const args = Array.prototype.slice.call(arguments);
		const msg = this.formatNoColor("SUC", args);
		this.append(msg);
		super.success(...args);
		return this;
	}
	
	warning ()
	{
		const args = Array.prototype.slice.call(arguments);
		const msg = this.formatNoColor("WAR", args);
		this.append(msg);
		super.warning(...args);
		return this;
	}
	
	warn ()
	{
		const args = Array.prototype.slice.call(arguments);
		const msg = this.formatNoColor("WAR", args);
		this.append(msg);
		super.warn(...args);
		return this;
	}
	
	error ()
	{
		const args = Array.prototype.slice.call(arguments);
		const msg = this.formatNoColor("ERR", args);
		this.append(msg);
		super.error(...args);
		return this;
	}
}

module.exports = p =>
{
	const l = new myLogger(require("path").resolve(p));
	
	const o_log = console.log;
	const o_warn = console.warn;
	const o_err = console.error;
	const o_debug = console.debug;
	const o_info = console.info;
	
	console.log = function () {
		const args = Array.prototype.slice.call(arguments);
		l.log(...args);
	}
	
	console.debug = function () {
		const args = Array.prototype.slice.call(arguments);
		l.debug(...args);
	}
	
	console.info = function () {
		const args = Array.prototype.slice.call(arguments);
		l.info(...args);
	}
	
	console.alert = function () {
		const args = Array.prototype.slice.call(arguments);
		l.alert(...args);
	}
	
	console.success = function () {
		const args = Array.prototype.slice.call(arguments);
		l.success(...args);
	}
	
	console.warning = function () {
		const args = Array.prototype.slice.call(arguments);
		l.warning(...args);
	}
	
	console.warn = function () {
		const args = Array.prototype.slice.call(arguments);
		l.warn(...args);
	}
	
	console.error = function () {
		const args = Array.prototype.slice.call(arguments);
		l.error(...args);
	}
	
	console.reset = function () {
		console.log = o_log;
		console.debug = o_debug;
		console.warn = o_warn;
		console.error = o_err;
		console.info = o_info;
		console.warning = undefined;
		console.success = undefined;
		console.alert = undefined;
	}
	
	return l;
}
