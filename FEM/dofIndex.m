function M = dofIndex(nodeDof, eleDof, ne, nc, method)

	switch method
	
		case 'contiguous'
			p = eleDof/nodeDof - 1;
			head = 1 + (0:ne-1)*(nodeDof*p);      % 1×(ne)
			if nc>1
				head = repelem(head, nc);         % 1×(ne*nc)
			end
			M = (0:eleDof-1).' + head;            % eleDof × (ne*nc)
	
	case 'strided'

		p = eleDof/nodeDof - 1;
		head = 1 + (0:ne-1)*(nodeDof*p);      % 1×(ne)

		M = (0:eleDof-1).' + head;            % eleDof × ne	
		if nc>1
			M = kron(ones(1,nc), M);      % eleDof × (ne*nc)
		end
	
	end
	
end
